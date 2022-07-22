#include "../src/cpp_macros.h"
!> \example simple_streamer.f90
!>
!> A simplified model for ionization waves and/or streamers in 2D
program streamer

  use m_af_types
  use m_af_core
  use m_af_ghostcell
  use m_af_utils
  use m_af_restrict
  use m_af_multigrid
  use m_af_output
  use m_write_silo
  use m_af_prolong

  implicit none

  integer            :: i, n
  character(len=100) :: fname
  type(af_t)         :: tree ! This contains the full grid information
  type(mg_t)         :: mg   ! Multigrid option struct
  type(mg_t)         :: ph_mg(3)
  type(ref_info_t)   :: refine_info


  ! Indices of cell-centered variables
  integer :: i_elec     ! Electron density
  integer :: i_pion     ! Positive ion density
  integer :: i_elec_old ! For time-stepping scheme
  integer :: i_pion_old ! For time-stepping scheme
  integer :: i_phi      ! Electrical potential
  integer :: i_fld      ! Electric field norm
  integer :: i_rhs      ! Source term Poisson

  !Cell centered variables related to photoionization
  integer :: i_ph_src !Photoionization source term
  integer :: i_ph_dist !Distribution of photoionizing photons

  integer :: i_psi_1 !Isotropic component of photon distribution function
  integer :: i_psi_2 !For the jth wavelength group
  integer :: i_psi_3 !j = 1, 2, 3

  integer :: i_rhs_1 !Right hand variable for Eddington Approximation
  integer :: i_rhs_2
  integer :: i_rhs_3

  integer :: i_tmp !Temporary variable

  ! Indices of face-centered variables **
  integer :: f_elec ! Electron flux
  integer :: f_fld  ! Electric field vector

  integer :: psi_indices(3) !arrays to hold indices of psi for easy access in loops
  integer :: rhs_indices(3) !same for rhs

  ! Simulation parameters
  real(dp), parameter :: end_time      = 10e-9_dp
  real(dp), parameter :: dt_output     = 20e-11_dp
  real(dp), parameter :: dt_max        = 20e-11_dp
  integer, parameter  :: ref_per_steps = 2
  integer, parameter  :: box_size      = 8

  ! Physical parameters
  real(dp), parameter :: applied_field = -1e7_dp
  real(dp), parameter :: mobility      = 0.03_dp
  real(dp), parameter :: diffusion_c   = 0.2_dp

  !photo-ionization parameters
  real(dp), parameter :: ph_quenching = 1
  real(dp), parameter :: ph_chi_min = 0.035 !/(torr*cm)
  real(dp), parameter :: ph_chi_max = 2 !/(torr*cm)
  real(dp), parameter :: I_0 = 3.5e22_dp !/(cm**3*s)
  real(dp), parameter :: p_o2 = 150 !Partial pressure of oxygen; 150 at ATP
  real(dp), parameter :: ui_ratio = 0.6 !Segur et. al. ss. 4.2.1
  real(dp), parameter :: p = 760 !torr; atmospheric pressure
  real(dp), parameter :: p_q = 30 !tore; see Segur. et. al. ss. 4.2.1
  real(dp), parameter :: xi = 0.2 !photoionization efficiency
  !Fit parameters
  !Reported in A Bourdon et al 2007 Plasma Sources Sci. Technol 16 656
  real(dp), parameter :: A(3) = (/ 0.0067_dp,0.0346_dp,0.3059_dp /)
  real(dp), parameter :: lambda(3) = (/ 0.0447_dp,0.1121_dp,0.5994_dp /)

  ! Computational domain
  real(dp), parameter :: domain_length = 2e-3_dp
  real(dp), parameter :: refine_max_dx = 1e-3_dp
  real(dp), parameter :: refine_min_dx = 1e-9_dp
  real(dp), parameter :: pi            = 3.14159265359
  ! Settings for the initial conditions
  real(dp), parameter :: current_density = 1e20_dp
  real(dp), parameter :: z_min   = 0.4_dp * domain_length
  real(dp), parameter :: z_max   = 0.6_dp * domain_length
  real(dp), parameter :: r_min = 0
  real(dp), parameter :: r_max = 0.1_dp * domain_length
  real(dp), parameter :: init_density = 1

  ! Simulation variables
  real(dp) :: dt
  real(dp) :: time
  integer  :: output_count

  call af_add_cc_variable(tree, "elec", ix=i_elec, n_copies=2)
  i_elec_old = af_find_cc_variable(tree, "elec_2")
  call af_add_cc_variable(tree, "pion", ix=i_pion, n_copies=2)
  i_pion_old = af_find_cc_variable(tree, "pion_2")
  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "fld", ix=i_fld)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "ph_src", ix=i_ph_src)
  call af_add_cc_variable(tree, "ph_dist", ix=i_ph_dist)
  call af_add_cc_variable(tree, "psi_1", ix=i_psi_1)
  call af_add_cc_variable(tree, "psi_2", ix=i_psi_2)
  call af_add_cc_variable(tree, "psi_3", ix=i_psi_3)
  call af_add_cc_variable(tree, "rhs_1", ix=i_rhs_1)
  call af_add_cc_variable(tree, "rhs_2", ix=i_rhs_2)
  call af_add_cc_variable(tree, "rhs_3", ix=i_rhs_3)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  call af_add_fc_variable(tree, "f_elec", ix=f_elec)
  call af_add_fc_variable(tree, "f_fld", ix=f_fld)

  !Initialize arrays for access to psi and rhs for each Eddington group

  psi_indices(1) = i_psi_1
  psi_indices(2) = i_psi_2
  psi_indices(3) = i_psi_3

  rhs_indices(1) = i_rhs_1
  rhs_indices(2) = i_rhs_2
  rhs_indices(3) = i_rhs_3

  ! Initialize the tree (which contains all the mesh information)
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [domain_length, domain_length], &
       [box_size, box_size], &
       coord=af_cyl, &
       periodic=[.false., .true.])

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs

  !Set multigrid options for photoionization PDEs
  ph_mg(1)%i_phi     = i_psi_1
  ph_mg(1)%i_rhs     = i_rhs_1
  ph_mg(1)%i_tmp     = i_tmp
  ph_mg(1)%helmholtz_lambda = (3 * (lambda(1) * p_o2) ** 2)

  ph_mg(2)%i_phi     = i_psi_2
  ph_mg(2)%i_rhs     = i_rhs_2
  ph_mg(2)%i_tmp     = i_tmp
  ph_mg(2)%helmholtz_lambda = (3 * (lambda(2) * p_o2) ** 2)

  ph_mg(3)%i_phi     = i_psi_3
  ph_mg(3)%i_rhs     = i_rhs_3
  ph_mg(3)%i_tmp     = i_tmp
  ph_mg(3)%helmholtz_lambda = (3 * (lambda(3) * p_o2) ** 2)

  ! Routines to use for ...
  mg%sides_bc => sides_bc_pot ! Filling ghost cell on physical boundaries

  !And for Helmholtz equations
  do i=1, 3
    if (i == 1) then
      ph_mg(i)%sides_bc => sides_bc
    else
      ph_mg(i)%sides_bc => af_bc_dirichlet_zero
    end if
    print *, "Helmholtz initialization"
    call mg_init(tree, ph_mg(i))
  end do

  ! This routine always needs to be called when using multigrid
  call mg_init(tree, mg)

  call af_set_cc_methods(tree, i_elec, af_bc_dirichlet_zero, &
       prolong=af_prolong_limit)
  call af_set_cc_methods(tree, i_pion, af_bc_dirichlet_zero, &
       prolong=af_prolong_limit)
  call af_set_cc_methods(tree, i_fld, af_bc_neumann_zero)

  output_count = 0 ! Number of output files written
  time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do i = 1, 4
    call af_adjust_refinement(tree, &
        init_refinement_routine, &
        refine_info)
     ! For each box in tree, set the initial conditions
     call af_loop_box(tree, set_initial_condition)

     ! Compute electric field on the tree.
     ! First perform multigrid to get electric potential,
     ! then take numerical gradient to geld field.
     call compute_fld(tree, .false.)

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     call af_adjust_refinement(tree, &               ! tree
          refinement_routine, & ! Refinement function
          refine_info)          ! Information about refinement

     ! If no new boxes have been added or removed, exit the loop
     if (refine_info%n_add == 0 .and. refine_info%n_rm == 0) exit
  end do

  call af_print_info(tree)

  do
     ! Get a new time step, which is at most dt_max
     call af_reduction(tree, &    ! Tree to do the reduction on
          max_dt, &  ! function
          get_min, & ! function
          dt_max, &  ! Initial value for the reduction
          dt)        ! Result of the reduction
      print *, "[*] Timestep", dt

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?", dt
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_count * dt_output <= time) then
        output_count = output_count + 1
        write(fname, "(A,I6.6)") "output/PEffStreamer_", output_count

        ! Write the cell centered data of a tree to a Silo file. Only the

        ! leaves of the tree are used
        call af_write_silo(tree, fname, output_count, time)
     end if

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ref_per_steps
        time = time + dt

        ! Copy previous solution
        call af_tree_copy_cc(tree, i_elec, i_elec_old)
        call af_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           ! First calculate fluxes
           call af_loop_tree(tree, fluxes_koren, leaves_only=.true.)
           call af_consistent_fluxes(tree, [f_elec])


           print *, "[*] Solving 3-Group Eddington for photoionization source term"
           !Solve 3-Group Eddington to get photoionization source
           call af_loop_box(tree, set_ph_dist)
           call compute_ph_src(tree, .true.)

           ! Update the solution
           call af_loop_box_arg(tree, update_solution, [dt], &
                leaves_only=.true.)

           ! Compute new field on first iteration
           if (i == 1) then
             call compute_fld(tree, .true.)
           end if
        end do

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        !call af_loop_box(tree, electron_beam)
        call af_loop_box(tree, average_dens)


        ! Compute field with new density
        call compute_fld(tree, .true.)
     end do

     ! Fill ghost cells for i_pion and i_elec
     call af_gc_tree(tree, [i_elec, i_pion])

     ! Adjust the refinement of a tree using refine_routine (see below) for grid
     ! refinement.
     ! Routine af_adjust_refinement sets the bit af_bit_new_children for each box
     ! that is refined.  On input, the tree should be balanced. On output,
     ! the tree is still balanced, and its refinement is updated (with at most
     ! one level per call).
     call af_adjust_refinement(tree, &               ! tree
          refinement_routine, & ! Refinement function
          refine_info, &        ! Information about refinement
          4)                    ! Buffer width (in cells)

     if (refine_info%n_add > 0 .or. refine_info%n_rm > 0) then
        ! Compute the field on the new mesh
        call compute_fld(tree, .true.)
     end if
  end do

contains

  subroutine init_refinement_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(box%n_cell, box%n_cell)
    integer                 :: i, j, nc
    real(dp)                :: rz(2)

    nc = box%n_cell

    do j=1, nc
      do i = 1, nc
        rz = af_r_cc(box, [i,j])

        !if (rz(2) > z_min .and. rz(2) < z_max) then
          !if (rz(1) > r_min .and. rz(1) < r_max) then
            !print *, "Refine"
            !cell_flags(i,j) = af_do_ref
          !else
            !print *, "Do not refine"
            !cell_flags(i,j) = af_keep_ref
        !  end if
        !end if
        cell_flags(i,j) = af_keep_ref
      end do
    end do
  end subroutine init_refinement_routine



  !> This routine sets the refinement flag for boxes(id)
  subroutine refinement_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    integer                  :: i, j, nc
    real(dp)                 :: dx, dens, fld, adx, rz(2)

    nc = box%n_cell
    dx = maxval(box%dr)

    do j = 1, nc
       do i = 1, nc
          rz   = af_r_cc(box, [i,j])
          dens = box%cc(i, j, i_elec)
          fld = box%cc(i, j, i_fld)
          adx = get_alpha(fld) * dx

          if (dens > 1e0_dp .and. adx > 0.8_dp) then
             cell_flags(i, j) = af_do_ref
          else if (dx < 1.25e-5_dp .and. adx < 0.1_dp) then
             cell_flags(i, j) = af_rm_ref
          else
             cell_flags(i, j) = af_keep_ref
          end if
       end do
    end do

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    else if (dx > 0.5_dp * refine_max_dx) then
       where(cell_flags == af_rm_ref) cell_flags = af_keep_ref
    end if

  end subroutine refinement_routine

  !> This routine sets the initial conditions for each box
  subroutine set_initial_condition(box)
    type(box_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: rz(2), normal_rands(2), vol

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          rz  = af_r_cc(box, [i,j])
          vol = (box%dr(1) * box%dr(2))**1.5_dp

          if (rz(2) > z_min .and. rz(2) < z_max) then
            if (rz(1) > r_min .and. rz(1) < r_max) then
             ! Approximate Poisson distribution with normal distribution
             normal_rands = two_normals(vol * init_density, &
                  sqrt(vol * init_density))
             ! Prevent negative numbers
             box%cc(i, j, i_elec) = abs(normal_rands(1)) / vol
          else
             box%cc(i, j, i_elec) = 0
           end if
          end if
       end do
    end do

    box%cc(:, :, i_pion) = box%cc(:, :, i_elec)
    box%cc(:, :, i_phi)  = 0 ! Inital potential set to zero

  end subroutine set_initial_condition

  !> Return two normal random variates
  !> http://en.wikipedia.org/wiki/Marsaglia_polar_method
  function two_normals(mean, sigma) result(rands)
    real(dp), intent(in) :: mean, sigma
    real(dp) :: rands(2), sum_sq

    do
       call random_number(rands)
       rands = 2 * rands - 1
       sum_sq = sum(rands**2)
       if (sum_sq < 1.0_dp .and. sum_sq > 0.0_dp) exit
    end do
    rands = rands * sqrt(-2 * log(sum_sq) / sum_sq)
    rands = mean + rands * sigma
  end function two_normals

  !> This function computes the minimum val of a and b
  real(dp) function get_min(a, b)
    real(dp), intent(in) :: a, b
    get_min = min(a, b)
  end function get_min

  !> Get maximum time step based on e.g. CFL criteria
  real(dp) function max_dt(box)
    type(box_t), intent(in) :: box
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    integer :: i, j, nc
    real(dp)               :: fld(2)
    real(dp)               :: dt_cfl, dt_dif, dt_drt

    nc = box%n_cell
    dt_cfl = dt_max

    do j = 1, nc
       do i = 1, nc
          fld(1) = 0.5_dp * (box%fc(i, j, 1, f_fld) + &
               box%fc(i+1, j, 1, f_fld))
          fld(2) = 0.5_dp * (box%fc(i, j, 2, f_fld) + &
               box%fc(i, j+1, 2, f_fld))

          ! The 0.5 is here because of the explicit trapezoidal rule
          dt_cfl = min(dt_cfl, 0.5_dp / sum(abs(fld * mobility) / box%dr))
       end do
    end do

    ! Dielectric relaxation time
    dt_drt = abs(UC_eps0 / (UC_elem_charge * mobility * &
         maxval(box%cc(1:nc, 1:nc, i_elec)) + epsilon(1.0_dp)))

    ! Diffusion condition
    dt_dif = 0.25_dp * minval(box%dr)**2 / max(diffusion_c, epsilon(1.0_dp))

    max_dt = min(dt_cfl, dt_drt, dt_dif)
    !print *, "cfl", dt_cfl, "drt", dt_drt, "dif", dt_dif
  end function max_dt


  !> This function gets the alpha value
  !>
  !> Taken from: Spatially hybrid computations for streamer discharges: II. Fully
  !> 3D simulations, Chao Li, Ute Ebert, Willem Hundsdorfer, J. Comput. Phys.
  !> 231, 1020-1050 (2012), doi:10.1016/j.jcp.2011.07.023
  elemental function get_alpha(fld) result(alpha)
    real(dp), intent(in) :: fld
    real(dp)             :: alpha
    real(dp), parameter  :: c0 = 1.04e1_dp
    real(dp), parameter  :: c1 = 0.601_dp
    real(dp), parameter  :: c2 = 1.86e7_dp

    alpha = exp(c0) * (abs(fld) * 1e-5_dp)**c1 * exp(-c2 / abs(fld))
  end function get_alpha

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, have_guess)
    type(af_t), intent(inout) :: tree
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = 1.6021766208e-19_dp / 8.8541878176e-12_dp
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(:, :, i_rhs) = fac * (&
               tree%boxes(id)%cc(:, :, i_elec) - &
               tree%boxes(id)%cc(:, :, i_pion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    ! Perform an FMG cycle (logicals: store residual, first call)
    call mg_fas_fmg(tree, mg, .false., have_guess)

    ! Compute field from potential
    call af_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call af_gc_tree(tree, [i_fld])
  end subroutine compute_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr(2)

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fc(1:nc+1, 1:nc, 1, f_fld) = inv_dr(1) * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, f_fld) = inv_dr(2) * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fc(1:nc, 1:nc, 1, f_fld) + &
         box%fc(2:nc+1, 1:nc, 1, f_fld))**2 + &
         0.25_dp * (box%fc(1:nc, 1:nc, 2, f_fld) + &
         box%fc(1:nc, 2:nc+1, 2, f_fld))**2)
  end subroutine fld_from_pot

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(tree, id)
    use m_af_flux_schemes
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: id
    real(dp)                  :: inv_dr(2)
    real(dp)                  :: cc(-1:tree%n_cell+2, -1:tree%n_cell+2, 1)
    real(dp), allocatable     :: v(:, :, :), dc(:, :, :)
    integer                   :: nc

    nc     = tree%n_cell
    inv_dr = 1/tree%boxes(id)%dr

    allocate(v(1:nc+1, 1:nc+1, 2))
    allocate(dc(1:nc+1, 1:nc+1, 2))

    call af_gc2_box(tree, id, [i_elec], cc)

    v = -mobility * tree%boxes(id)%fc(:, :, :, f_fld)
    dc = diffusion_c

    call flux_koren_2d(cc(:, :, 1), v, nc, 2)
    call flux_diff_2d(cc(:, :, 1), dc, inv_dr, nc, 2)

    tree%boxes(id)%fc(:, :, :, f_elec) = v + dc
  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr(2), src, sflux, fld, ph_src
    real(dp)                    :: alpha
    real(dp)                    :: beam_src
    integer                     :: i, j, nc

    nc                    = box%n_cell
    inv_dr                = 1/box%dr

    do j = 1, nc
       do i = 1, nc
          fld   = box%cc(i,j, i_fld)
          alpha = get_alpha(fld)
          beam_src = get_beam_src(box,i,j)
          src   = box%cc(i, j, i_elec) * mobility * abs(fld) * alpha
          ph_src = box%cc(i, j, i_ph_src)

          sflux = inv_dr(1) * (box%fc(i, j, 1, f_elec) - &
               box%fc(i+1, j, 1, f_elec)) + &
               inv_dr(2) * (box%fc(i, j, 2, f_elec) - &
               box%fc(i, j+1, 2, f_elec))

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (beam_src + src + sflux + ph_src) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + (src + ph_src) * dt(1)
       end do
    end do

  end subroutine update_solution

  !Return the source due to the electron beam
  function get_beam_src(box,i,j) result(beam_src)
    integer, intent(in)     :: i,j
    type(box_t), intent(in) :: box
    real(dp)                :: beam_src
    real(dp)                :: rz(2), dr(2), vol, normal_rands(2)

    rz = af_r_cc(box, [i,j])
    dr = box%dr
    vol = (dr(1) * dr(2))**1.5

    !If the beam source is within the beam source region, add density there
    !Randomize that density a bit to be more realistic
    if (rz(2) > z_min .and. rz(2) < z_max) then
      if (rz(1) > r_min .and. rz(1) < r_max) then
        normal_rands = two_normals(vol * current_density, &
           sqrt(vol * current_density))
           ! Prevent negative numbers
           beam_src = abs(normal_rands(1)) / vol
         else
           beam_src = 0
         end if
       end if
  end function get_beam_src

  !> This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type

    select case (nb)
    case (af_neighb_lowx)
       bc_type = af_bc_neumann
       bc_val  = 0.0_dp
    case (af_neighb_highx)
       bc_type = af_bc_dirichlet
       bc_val = applied_field
    case default
       stop "sides_bc_pot: unspecified boundary"
    end select
  end subroutine sides_bc_pot

  !Subroutines for photoionization
  subroutine compute_ph_src(tree, have_guess)
    type(af_t), intent(inout)   :: tree
    logical, intent(in)         :: have_guess
    real(dp), parameter         :: rhs_c = -3 * p_o2 / (xi * 2.99792458e8) !Constant on right hand side of Helmholtz equations; -3p_o2/c
    integer                     :: lvl, i, j, k, id, nc

    nc = tree%n_cell

    do j = 1, 3
      !set the source term (rhs)
      do lvl = 1, tree%highest_lvl
        do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
            tree%boxes(id)%cc(:, :, rhs_indices(j)) = rhs_c * lambda(j) * (&
                 tree%boxes(id)%cc(:, :, i_ph_dist))
        end do
      end do
      !Perform an FMG cycle for each j=1,2,3
      call mg_fas_fmg(tree, ph_mg(j), .false., have_guess)
    end do

    call af_loop_box(tree, ph_src_from_psi)
  end subroutine compute_ph_src

  subroutine ph_src_from_psi(box)
    type(box_t), intent(inout)  :: box
    integer                     :: nc, i, j, k
    real(dp)                    :: s_ph

    nc = box%n_cell

    do i = 1, nc
      do j = 1, nc
        s_ph = 0_dp

        !Update s_ph using equation (20) from Bourdon et al
        do k = 1, 3
          s_ph = s_ph + A(k) * p_o2 * 2.99792458e8 * box%cc(i, j, psi_indices(k))
        end do

        box%cc(i, j, i_ph_src) = s_ph
    end do
  end do
  end subroutine ph_src_from_psi

  !Set the boundary conditions for the smallest lambda
  !Using the Zhelensky model
  subroutine sides_bc(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(2, box%n_cell)
    real(dp), intent(out)   :: bc_val(box%n_cell)
    integer, intent(out)    :: bc_type
    integer                 :: nc, lvl, i, j, k, l, m, id
    real(dp)                :: s_ph, rz_1(2), rz_2(2), R, vol, dr(2)

    !Dirichlet boundary conditions
    bc_type = af_bc_dirichlet


    do j = 1, box%n_cell
      rz_1(1) = coords(1, j)
      rz_1(2) = coords(2, j)
      s_ph = 0

      do lvl = 1, tree%highest_lvl
        do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
          !print *, "Level, i: ", lvl, i

          !Solve for boundary conditions using Zheleznyak integral
          do l = 1, nc
            do m = 1, nc

              rz_2 = af_r_cc(tree%boxes(id), [l, m])
              R = sqrt((rz_2(1) - rz_1(1))**2 + (rz_2(2) - rz_1(2))**2)
              !print *, R, rz_1
              dr = tree%boxes(id)%dr
              vol = (dr(1) * dr(2)) ** 1.5

              s_ph = s_ph + tree%boxes(id)%cc(l, m, i_ph_dist) * vol * &
                        f(R) / (4 * pi * R**2)
            end do
          end do

        end do
      end do

      !print *, s_ph
      bc_val(j) = s_ph / (p_o2 * A(1) * 2.99792458e8)
    end do
    !print *, bc_val
  end subroutine sides_bc

  subroutine set_ph_dist(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc, i, j

    nc = box%n_cell
    do i = 1, nc
      do j = 1, nc
        box%cc(i, j, i_ph_dist) = p_q / (p + p_q) * ui_ratio * box%cc(i, j, i_elec)
      end do
    end do
  end subroutine set_ph_dist

  real(dp) function f(R)
    real(dp), intent(in) :: R

    f = (exp(-ph_chi_min * p_o2 * R) - exp(-ph_chi_max * p_o2 * R)) / (R * log(ph_chi_max / ph_chi_min))
  end function f
end program streamer
