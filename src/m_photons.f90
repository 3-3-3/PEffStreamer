!> Module that provides routines for operations on photons
module m_photons

  use m_lookup_table

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  type photoi_tbl_t
     type(lookup_table_t) :: tbl           !< The lookup table
     real(dp)             :: frac_in_tbl   !< Fraction photons in table
  end type photoi_tbl_t

  public :: photoi_tbl_t

  ! Public methods
  public :: photoi_print_grid_spacing
  public :: photoi_get_table_air
  public :: photoi_do_absorption
  public :: photoi_absorption_func_air
  public :: photoi_set_src_2d
  public :: photoi_set_src_3d

contains

  !> Print the grid spacing used for the absorption of photons
  subroutine photoi_print_grid_spacing(pi_tbl, dr_base, fac_dx)
    type(photoi_tbl_t), intent(in) :: pi_tbl
    real(dp), intent(in)           :: dr_base
    real(dp), intent(in)           :: fac_dx
    real(dp)                       :: lengthscale, dx
    integer                        :: lvl

    ! Get a typical length scale for the absorption of photons
    lengthscale = LT_get_col(pi_tbl%tbl, 1, fac_dx)

    ! Determine at which level we estimate the photoionization source term. This
    ! depends on the typical length scale for absorption.
    lvl = get_lvl_length(dr_base, lengthscale)
    dx  = dr_base * (0.5_dp**(lvl-1))

    write(*, "(A,E12.3)") " Photoionization spacing: ", dx
  end subroutine photoi_print_grid_spacing

  !> Compute the photonization table for air. If the absorption function is
  !> f(r), this table contains r as a function of F (the cumulative absorption
  !> function, between 0 and 1). Later on a distance can be sampled by drawing a
  !> uniform 0,1 random number and finding the corresponding distance. The table
  !> constructed up to max_dist; we can ignore photons that fly very far.
  subroutine photoi_get_table_air(photoi_tbl, p_O2, max_dist)
    !< The photonization table
    type(photoi_tbl_t), intent(inout) :: photoi_tbl
    !< Partial pressure of oxygen (bar)
    real(dp), intent(in)              :: p_O2
    !< Maximum distance in lookup table
    real(dp), intent(in)  :: max_dist

    integer               :: n
    integer, parameter    :: tbl_size  = 500
    real(dp), allocatable :: fsum(:), dist(:)
    real(dp)              :: dF, drdF, r, F, Fmax_guess

    ! First estimate which fraction of photons are within max_dist (start with
    ! the upper bound of 1.0)
    Fmax_guess = 1.0_dp

    ! 5 loops should be enough for a good guess
    do n = 1, 5
       ! When Fmax_guess decreases so will dF, since we keep the number of
       ! points the same
       dF = Fmax_guess / (tbl_size-1)
       r  = 0
       F  = 0

       do
          drdF = rk4_drdF(r, dF, p_O2)
          r = r + dF * drdF
          F = F + df

          if (r > max_dist) then
             ! A better estimate for the upper bound
             Fmax_guess = F
             exit
          end if
       end do
    end do

    ! Make arrays larger so that we surely are in range of maxdist
    allocate(fsum(2 * tbl_size))
    allocate(dist(2 * tbl_size))

    ! Now create table
    dF = Fmax_guess / (tbl_size-1)
    dist(1) = 0
    fsum(1) = 0

    ! Compute r(F) for F = dF, 2 dF, 3 dF, ...
    do n = 2, 2 * tbl_size
       drdF = rk4_drdF(dist(n-1), dF, p_O2)
       fsum(n) = fsum(n-1) + dF
       dist(n) = dist(n-1) + dF * drdF
       if (dist(n) > max_dist) exit
    end do

    if (n > tbl_size + 10) &
         stop "photoi_get_table_air: integration accuracy fail"

    ! Scale table to lie between 0 and 1 (later on we have to correct for this)
    photoi_tbl%frac_in_tbl = fsum(n-1)
    fsum(1:n-1) = fsum(1:n-1) / fsum(n-1)

    photoi_tbl%tbl = LT_create(0.0_dp, 1.0_dp, tbl_size, 1)
    call LT_set_col(photoi_tbl%tbl, 1, fsum(1:n-1), dist(1:n-1))
  end subroutine photoi_get_table_air

  !> Runge Kutta 4 method
  real(dp) function rk4_drdF(r, dF, p_O2)
    real(dp), intent(in) :: r    !> Initial point
    real(dp), intent(in) :: dF   !> grid size
    real(dp), intent(in) :: p_O2 !< Partial pressure of oxygen (bar)
    real(dp)             :: drdF
    real(dp)             :: sum_drdF
    real(dp), parameter  :: one_sixth = 1 / 6.0_dp

    ! Step 1 (at initial r)
    drdF = 1 / photoi_absorption_func_air(r, p_O2)
    sum_drdF = drdF

    ! Step 2 (at initial r + dr/2)
    drdF = 1 / photoi_absorption_func_air(r + 0.5_dp * dF * drdF, p_O2)
    sum_drdF = sum_drdF + 2 * drdF

    ! Step 3 (at initial r + dr/2)
    drdF = 1 / photoi_absorption_func_air(r + 0.5_dp * dF * drdF, p_O2)
    sum_drdF = sum_drdF + 2 * drdF

    ! Step 4 (at initial r + dr)
    drdF = 1 / photoi_absorption_func_air(r + dF * drdF, p_O2)
    sum_drdF = sum_drdF + drdF

    ! Combine r derivatives at steps
    rk4_drdF = one_sixth * sum_drdF
  end function rk4_drdF

  !> The absorption function for photons in air according to Zheleznyak's model
  real(dp) function photoi_absorption_func_air(dist, p_O2)
    use m_units_constants
    real(dp), intent(in) :: dist   !< Distance
    real(dp), intent(in) :: p_O2   !< Partial pressure of oxygen (bar)
    real(dp)             :: r
    real(dp), parameter  :: c0 = 3.5_dp / UC_torr_to_bar
    real(dp), parameter  :: c1 = 200 / UC_torr_to_bar
    real(dp), parameter  :: eps = epsilon(1.0_dp)

    r = p_O2 * dist
    if (r * (c0 + c1) < eps) then
       ! Use limit to prevent over/underflow
       photoi_absorption_func_air = (c1 - c0 + 0.5_dp * (c0**2 - c1**2) * r) &
            * p_O2 / log(c1/c0)
    else if (r * c0 > -log(eps)) then
       ! Use limit to prevent over/underflow
       photoi_absorption_func_air = eps
    else
       photoi_absorption_func_air = (exp(-c0 * r) - exp(-c1 * r)) / (dist * log(c1/c0))
    end if
  end function photoi_absorption_func_air

  !> Determine the lowest level at which the grid spacing is smaller than 'length'.
  integer function get_lvl_length(dr_base, length)
    real(dp), intent(in) :: dr_base !< cell spacing at lvl 1
    real(dp), intent(in) :: length  !< Some length
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    real(dp) :: ratio

    ratio = dr_base / length
    if (ratio <= 1) then
       get_lvl_length = 1
    else
       get_lvl_length = 1 + ceiling(log(ratio) * invlog2)
    end if
  end function get_lvl_length

  !> As get_lvl_length but with a random choice between lvl and lvl-1
  integer function get_rlvl_length(dr_base, length, rng)
    use m_random
    real(dp), intent(in) :: dr_base   !< cell spacing at lvl 1
    real(dp), intent(in) :: length    !< Some length
    type(RNG_t), intent(inout) :: rng !< Random number generator
    real(dp), parameter :: invlog2 = 1 / log(2.0_dp)
    real(dp) :: ratio, tmp

    ratio = dr_base / length
    if (ratio <= 1) then
       get_rlvl_length = 1
    else
       tmp = log(ratio) * invlog2
       get_rlvl_length = floor(tmp)
       if (rng%unif_01() < tmp - get_rlvl_length) &
            get_rlvl_length = get_rlvl_length + 1
    end if
  end function get_rlvl_length

  !> Given a list of photon production positions (xyz_in), compute where they
  !> end up (xyz_out).
  subroutine photoi_do_absorption(xyz_in, xyz_out, n_dim, n_photons, tbl, rng)
    use m_lookup_table
    use m_random
    use omp_lib
    integer, intent(in)              :: n_photons !< Number of photons
    !< Input (x,y,z) values
    real(dp), intent(in)             :: xyz_in(3, n_photons)
    !< Output (x,y,z) values
    real(dp), intent(out)            :: xyz_out(3, n_photons)
    integer, intent(in)              :: n_dim     !< 2 or 3 dimensional
    !< Lookup table
    type(lookup_table_t), intent(in) :: tbl
    type(RNG_t), intent(inout)       :: rng       !< Random number geneator
    integer                          :: n, n_procs, proc_id
    real(dp)                         :: rr, dist
    type(PRNG_t)                     :: prng

    !$omp parallel private(n, rr, dist, proc_id)
    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init_parallel(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()

    if (n_dim == 2) then
       !$omp do
       do n = 1, n_photons
          rr = prng%rngs(proc_id)%unif_01()
          dist = LT_get_col(tbl, 1, rr)
          xyz_out(1:n_dim, n) =  xyz_in(1:n_dim, n) + &
               prng%rngs(proc_id)%circle(dist)
       end do
       !$omp end do
    else if (n_dim == 3) then
       !$omp do
       do n = 1, n_photons
          rr = prng%rngs(proc_id)%unif_01()
          dist = LT_get_col(tbl, 1, rr)
          xyz_out(:, n) =  xyz_in(:, n) + prng%rngs(proc_id)%sphere(dist)
       end do
       !$omp end do
    else
       print *, "photoi_do_absorption: unknown n_dim", n_dim
       stop
    end if
    !$omp end parallel
  end subroutine photoi_do_absorption

  !> Set the source term due to photoionization for 2D models. At most
  !> max_photons discrete photons are produced.
  subroutine photoi_set_src_2d(tree, pi_tbl, rng, max_photons, &
       i_src, i_photo, fac_dx, const_dx, use_cyl, min_dx, dt)
    use m_random
    use m_a2_types
    use m_a2_utils
    use m_a2_ghostcell
    use m_a2_prolong
    use m_lookup_table
    use omp_lib

    type(a2_t), intent(inout)  :: tree   !< Tree
    type(photoi_tbl_t)         :: pi_tbl !< Table to sample abs. lengths
    type(RNG_t), intent(inout) :: rng    !< Random number generator
    !> Maximum number of discrete photons to use
    integer, intent(in)        :: max_photons
    !> Input variable that contains photon production per cell
    integer, intent(in)        :: i_src
    !> Input variable that contains photoionization source rate
    integer, intent(in)        :: i_photo
    !> Use dx proportional to this value
    real(dp), intent(in)       :: fac_dx
    !> Use constant grid spacing or variable
    logical, intent(in)        :: const_dx
    !> Use cylindrical coordinates
    logical, intent(in)        :: use_cyl
    !> Minimum spacing for absorption
    real(dp), intent(in)        :: min_dx
    !> Time step, if present use "physical" photons
    real(dp), intent(in), optional :: dt

    integer                     :: lvl, ix, id, nc, min_lvl, highest_lvl
    integer                     :: i, j, n, n_create, n_used, i_ph
    integer                     :: proc_id, n_procs
    integer                     :: pho_lvl
    real(dp)                    :: tmp, dr, fac, dist, r(3)
    real(dp)                    :: sum_production, pi_lengthscale
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_abs(:, :)
    real(dp), parameter         :: pi = acos(-1.0_dp)
    type(PRNG_t)                :: prng
    type(a2_loc_t), allocatable :: ph_loc(:)

    nc = tree%n_cell

    ! Compute the sum of photon production
    call a2_tree_sum_cc(tree, i_src, sum_production)

    if (present(dt)) then
       ! Create "physical" photons when less than max_photons are produced
       fac = min(dt, max_photons / (sum_production + epsilon(1.0_dp)))
    else
       ! Create approximately max_photons
       fac = max_photons / (sum_production + epsilon(1.0_dp))
    end if

    ! Allocate a bit more space because of stochastic production
    allocate(xyz_src(3, nint(1.2_dp * fac * sum_production + 1000)))

    ! Now loop over all leaves and create photons using random numbers
    n_used = 0

    !$omp parallel private(lvl, ix, id, i, j, n, r, dr, i_ph, proc_id, &
    !$omp tmp, n_create)

    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init_parallel(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()

    do lvl = 1, tree%highest_lvl
       dr = a2_lvl_dr(tree, lvl)
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do j = 1, nc
             do i = 1, nc
                if (tree%boxes(id)%coord_t == af_cyl) then
                   tmp = a2_cyl_radius_cc(tree%boxes(id), i)
                   tmp = fac * 2 * pi * tmp * &
                        tree%boxes(id)%cc(i, j, i_src) * dr**2
                else
                   tmp = fac * tree%boxes(id)%cc(i, j, i_src) * dr**2
                end if

                n_create = floor(tmp)

                if (prng%rngs(proc_id)%unif_01() < tmp - n_create) &
                     n_create = n_create + 1

                if (n_create > 0) then
                   !$omp critical
                   i_ph = n_used
                   n_used = n_used + n_create
                   !$omp end critical

                   ! Location of production
                   ! TODO: not all from cell center
                   r(1:2) = a2_rr_cc(tree%boxes(id), [i-0.5_dp, j+0.0_dp])
                   r(3) = 0

                   do n = 1, n_create
                      xyz_src(:, i_ph+n) = r
                   end do
                end if
             end do
          end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    allocate(xyz_abs(3, n_used))
    allocate(ph_loc(n_used))


    if (use_cyl) then
       ! Get location of absorbption
       call photoi_do_absorption(xyz_src, xyz_abs, 3, n_used, pi_tbl%tbl, rng)

       !$omp do
       do n = 1, n_used
          ! Set x coordinate to radius (norm of 1st and 3rd coord.)
          xyz_abs(1, n) = sqrt(xyz_abs(1, n)**2 + xyz_abs(3, n)**2)
       end do
       !$omp end do
    else
       ! Get location of absorbption
       call photoi_do_absorption(xyz_src, xyz_abs, 2, n_used, pi_tbl%tbl, rng)
    end if

    if (const_dx) then
       ! Get a typical length scale for the absorption of photons
       pi_lengthscale = LT_get_col(pi_tbl%tbl, 1, fac_dx)

       ! Determine at which level we estimate the photoionization source term. This
       ! depends on the typical length scale for absorption.
       pho_lvl = get_lvl_length(tree%dr_base, pi_lengthscale)

       !$omp parallel do
       do n = 1, n_used
          ph_loc(n) = a2_get_loc(tree, xyz_abs(1:2, n), pho_lvl)
       end do
       !$omp end parallel do
    else
       highest_lvl = get_lvl_length(tree%dr_base, min_dx)
       !$omp parallel private(n, dist, lvl, proc_id)
       proc_id = 1+omp_get_thread_num()
       !$omp do
       do n = 1, n_used
          dist = norm2(xyz_abs(1:2, n) - xyz_src(1:2, n))
          lvl = get_rlvl_length(tree%dr_base, fac_dx * dist, prng%rngs(proc_id))
          if (lvl > highest_lvl) lvl = highest_lvl
          ph_loc(n) = a2_get_loc(tree, xyz_abs(1:2, n), lvl)
       end do
       !$omp end do
       !$omp end parallel
    end if

    ! Clear variable i_photo, in which we will store the photoionization source term

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a2_box_clear_cc(tree%boxes(id), i_photo)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Add photons to production rate. Currently, this is done sequentially.
    if (use_cyl) then
       tmp = fac * 2 * pi       ! Temporary factor to speed up loop

       do n = 1, n_used
          id = ph_loc(n)%id
          if (id > af_no_box) then
             i = ph_loc(n)%ix(1)
             j = ph_loc(n)%ix(2)
             dr = tree%boxes(id)%dr
             r(1:2) = a2_r_cc(tree%boxes(id), [i, j])
             tree%boxes(id)%cc(i, j, i_photo) = &
                  tree%boxes(id)%cc(i, j, i_photo) + &
                  pi_tbl%frac_in_tbl/(tmp * dr**2 * r(1))
          end if
       end do
    else
       do n = 1, n_used
          id = ph_loc(n)%id
          if (id > af_no_box) then
             i = ph_loc(n)%ix(1)
             j = ph_loc(n)%ix(2)
             dr = tree%boxes(id)%dr
             tree%boxes(id)%cc(i, j, i_photo) = &
                  tree%boxes(id)%cc(i, j, i_photo) + &
                  pi_tbl%frac_in_tbl/(fac * dr**2)
          end if
       end do
    end if

    ! Set ghost cells on highest level with photon source
    if (const_dx) then
       min_lvl = pho_lvl
    else
       min_lvl = 1
    end if

    !$omp parallel private(lvl, i, id)
    ! Prolong to finer grids
    do lvl = min_lvl, tree%highest_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_gc_box(tree%boxes, id, i_photo, &
               a2_gc_interp, a2_bc_neumann_zero)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a2_prolong_linear_from(tree%boxes, id, i_photo, add=.true.)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine photoi_set_src_2d

  subroutine photoi_set_src_3d(tree, pi_tbl, rng, max_photons, &
       i_src, i_photo, fac_dx, const_dx, min_dx, dt)
    use m_random
    use m_a3_types
    use m_a3_utils
    use m_a3_ghostcell
    use m_a3_prolong
    use m_lookup_table
    use omp_lib

    type(a3_t), intent(inout)   :: tree   !< Tree
    type(photoi_tbl_t)          :: pi_tbl !< Table to sample abs. lengths
    type(RNG_t), intent(inout)  :: rng    !< Random number generator
    !> Maximum number of discrete photons to use
    integer, intent(in)         :: max_photons
    !> Input variable that contains photon production per cell
    integer, intent(in)         :: i_src
    !> Input variable that contains photoionization source rate
    integer, intent(in)         :: i_photo
    !> Use dx proportional to this value
    real(dp), intent(in)        :: fac_dx
    !> Use constant grid spacing or variable
    logical, intent(in)         :: const_dx
    !> Minimum spacing for absorption
    real(dp), intent(in)        :: min_dx
    !> Time step, if present use "physical" photons
    real(dp), intent(in), optional :: dt

    integer                     :: lvl, ix, id, nc
    integer                     :: i, j, k, n, n_create, n_used, i_ph
    integer                     :: proc_id, n_procs
    integer                     :: pho_lvl, highest_lvl, min_lvl
    real(dp)                    :: tmp, dr, fac, dist
    real(dp)                    :: sum_production, pi_lengthscale
    real(dp), allocatable       :: xyz_src(:, :)
    real(dp), allocatable       :: xyz_abs(:, :)
    type(PRNG_t)                :: prng
    type(a3_loc_t), allocatable :: ph_loc(:)

    nc = tree%n_cell

    ! Compute the sum of photon production
    call a3_tree_sum_cc(tree, i_src, sum_production)

    if (present(dt)) then
       ! Create "physical" photons when less than max_photons are produced
       fac = min(dt, max_photons / (sum_production + epsilon(1.0_dp)))
    else
       ! Create approximately max_photons
       fac = max_photons / (sum_production + epsilon(1.0_dp))
    end if

    ! Allocate a bit more space because of stochastic production
    allocate(xyz_src(3, nint(1.2_dp * fac * sum_production + 1000)))

    ! Now loop over all leaves and create photons using random numbers
    n_used = 0

    !$omp parallel private(lvl, ix, id, i, j, k, n, dr, i_ph, &
    !$omp proc_id, tmp, n_create)
    !$omp single
    n_procs = omp_get_num_threads()
    call prng%init_parallel(n_procs, rng)
    !$omp end single

    proc_id = 1+omp_get_thread_num()
    tmp = 0
    do lvl = 1, tree%highest_lvl
       dr = a3_lvl_dr(tree, lvl)
       !$omp do
       do ix = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(ix)

          do k = 1, nc
             do j = 1, nc
                do i = 1, nc
                   tmp = fac * tree%boxes(id)%cc(i, j, k, i_src) * dr**3
                   n_create = floor(tmp)

                   if (prng%rngs(proc_id)%unif_01() < tmp - n_create) &
                        n_create = n_create + 1

                   if (n_create > 0) then
                      !$omp critical
                      i_ph = n_used
                      n_used = n_used + n_create
                      !$omp end critical

                      do n = 1, n_create
                         xyz_src(:, i_ph+n) = &
                              a3_r_cc(tree%boxes(id), [i, j, k])
                      end do
                   end if
                end do
             end do
          end do
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    allocate(xyz_abs(3, n_used))
    allocate(ph_loc(n_used))

    ! Get location of absorption
    call photoi_do_absorption(xyz_src, xyz_abs, 3, n_used, pi_tbl%tbl, rng)

    if (const_dx) then
       ! Get a typical length scale for the absorption of photons
       pi_lengthscale = LT_get_col(pi_tbl%tbl, 1, fac_dx)

       ! Determine at which level we estimate the photoionization source term. This
       ! depends on the typical length scale for absorption.
       pho_lvl = get_lvl_length(tree%dr_base, pi_lengthscale)

       !$omp parallel do
       do n = 1, n_used
          ph_loc(n) = a3_get_loc(tree, xyz_abs(:, n), pho_lvl)
       end do
       !$omp end parallel do
    else
       highest_lvl = get_lvl_length(tree%dr_base, min_dx)
       !$omp parallel private(n, dist, lvl, proc_id)
       proc_id = 1+omp_get_thread_num()
       !$omp do
       do n = 1, n_used
          dist = norm2(xyz_abs(:, n) - xyz_src(:, n))
          lvl = get_rlvl_length(tree%dr_base, fac_dx * dist, prng%rngs(proc_id))
          if (lvl > highest_lvl) lvl = highest_lvl
          ph_loc(n) = a3_get_loc(tree, xyz_abs(:, n), lvl)
       end do
       !$omp end do
       !$omp end parallel
    end if

    ! Clear variable i_photo, in which we will store the photoionization source term

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)
          call a3_box_clear_cc(tree%boxes(id), i_photo)
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Add photons to production rate. Currently, this is done sequentially.
    do n = 1, n_used
       id = ph_loc(n)%id
       if (id > af_no_box) then
          i = ph_loc(n)%ix(1)
          j = ph_loc(n)%ix(2)
          k = ph_loc(n)%ix(3)
          dr = tree%boxes(id)%dr
          tree%boxes(id)%cc(i, j, k, i_photo) = &
               tree%boxes(id)%cc(i, j, k, i_photo) + &
               pi_tbl%frac_in_tbl/(fac * dr**3)
       end if
    end do

    ! Set ghost cells on highest level with photon source
    if (const_dx) then
       min_lvl = pho_lvl
    else
       min_lvl = 1
    end if

    !$omp parallel private(lvl, i, id)
    ! Prolong to finer grids
    do lvl = min_lvl, tree%highest_lvl-1
       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a3_gc_box(tree%boxes, id, i_photo, &
               a3_gc_interp, a3_bc_neumann_zero)
       end do
       !$omp end do

       !$omp do
       do i = 1, size(tree%lvls(lvl)%parents)
          id = tree%lvls(lvl)%parents(i)
          call a3_prolong_linear_from(tree%boxes, id, i_photo, add=.true.)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine photoi_set_src_3d

end module m_photons
