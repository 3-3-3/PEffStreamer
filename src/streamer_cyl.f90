program streamer_cyl
  use m_afivo_2d
  use m_mg_2d
  use m_write_silo
  use m_lookup_table
  use m_config
  use m_random
  use m_photons

  implicit none

  integer, parameter :: dp       = kind(0.0d0)
  integer, parameter :: name_len = 200

  ! Indices of cell-centered variables
  integer, parameter :: n_var_cell = 9
  integer, parameter :: i_elec     = 1 ! Electron density
  integer, parameter :: i_pion     = 2 ! Positive ion density
  integer, parameter :: i_elec_old = 3 ! For time-stepping scheme
  integer, parameter :: i_pion_old = 4 ! For time-stepping scheme
  integer, parameter :: i_phi      = 5 ! Electrical potential
  integer, parameter :: i_fld      = 6 ! Electric field norm
  integer, parameter :: i_rhs      = 7 ! Source term Poisson
  integer, parameter :: i_pho      = 8 ! Phototionization rate
  integer, parameter :: i_eps      = 9 ! Level set function
  character(len=10)  :: cc_names(n_var_cell) = &
       [character(len=10) :: "elec", "pion", "elec_old", &
       "pion_old", "phi", "fld", "rhs", "pho", "eps"]

  ! Indices of face-centered variables
  integer, parameter :: n_var_face = 2
  integer, parameter :: f_elec     = 1 ! Electron flux
  integer, parameter :: f_fld      = 2 ! Electric field vector

  ! Indices of transport data
  integer, parameter :: n_var_td = 4
  integer, parameter :: i_mobility = 1
  integer, parameter :: i_diffusion = 2
  integer, parameter :: i_alpha = 3
  integer, parameter :: i_eta = 4

  character(len=name_len) :: sim_name, output_dir
  character(len=name_len) :: cfg_name, tmp_name, prev_name

  type initcnd_t
     real(dp)              :: bg_dens
     integer               :: n_cond
     real(dp), allocatable :: seed_r0(:, :)
     real(dp), allocatable :: seed_r1(:, :)
     real(dp), allocatable :: seed_dens(:)
     real(dp), allocatable :: seed_width(:)
     integer, allocatable  :: seed_falloff(:)
  end type initcnd_t

  type(initcnd_t)   :: init_cond
  type(LT_table_t)  :: td_tbl             ! Table with transport data vs fld
  type(CFG_t)       :: sim_cfg            ! The configuration for the simulation
  type(a2_t)        :: tree               ! This contains the full grid information
  type(mg2_t)       :: mg                 ! Multigrid option struct
  type(RNG_t)       :: sim_rng            ! Random number generator
  type(ref_info_t)  :: ref_info

  logical           :: photoi_enabled     ! Whether we use phototionization
  real(dp)          :: photoi_frac_O2     ! Oxygen fraction
  real(dp)          :: photoi_eta         ! Photoionization efficiency
  integer           :: photoi_num_photons ! Number of photons to use
  type(PH_tbl_t)    :: photoi_tbl         ! Table for photoionization

  real(dp) :: fld_mod_t0, fld_sin_amplitude, fld_sin_freq, fld_lin_deriv

  integer           :: i, n, n_steps_amr
  integer           :: output_cnt
  real(dp)          :: dt, time, end_time
  real(dp)          :: dt_output, dt_max
  character(len=200) :: fname, fname_axis, fname_stats
  logical           :: write_out

  ! How many multigrid FMG cycles we perform per time step
  integer, parameter :: n_fmg_cycles = 1

  ! The size of the boxes that we use to construct our mesh
  integer :: box_size

  ! The length of the (square) domain
  real(dp) :: domain_len

  ! The applied electric field
  real(dp) :: applied_fld

  ! The applied voltage
  real(dp) :: applied_voltage

  ! Pressure of the gas in bar
  real(dp) :: gas_pressure

  ! Dielectric constant
  real(dp) :: epsilon_diel

  call create_cfg(sim_cfg)

  sim_name = ""
  prev_name = ""
  do n = 1, command_argument_count()
     call get_command_argument(n, cfg_name)
     call CFG_read_file(sim_cfg, trim(cfg_name))

     call CFG_get(sim_cfg, "sim_name", tmp_name)
     if (sim_name == "") then
        sim_name = tmp_name
     else if (tmp_name /= "" .and. tmp_name /= prev_name) then
        sim_name = trim(sim_name) // "_" // trim(tmp_name)
     end if
     prev_name = tmp_name
  end do

  call CFG_get(sim_cfg, "end_time", end_time)
  call CFG_get(sim_cfg, "box_size", box_size)
  call CFG_get(sim_cfg, "output_dir", output_dir)
  call CFG_get(sim_cfg, "domain_len", domain_len)
  call CFG_get(sim_cfg, "applied_fld", applied_fld)
  call CFG_get(sim_cfg, "fld_mod_t0", fld_mod_t0)
  call CFG_get(sim_cfg, "fld_sin_amplitude", fld_sin_amplitude)
  call CFG_get(sim_cfg, "fld_sin_freq", fld_sin_freq)
  call CFG_get(sim_cfg, "fld_lin_deriv", fld_lin_deriv)
  call CFG_get(sim_cfg, "dt_output", dt_output)
  call CFG_get(sim_cfg, "num_steps_amr", n_steps_amr)
  call CFG_get(sim_cfg, "dt_max", dt_max)
  call CFG_get(sim_cfg, "epsilon_diel", epsilon_diel)

  if (trim(output_dir) == "") stop "No output directory given"
  tmp_name = trim(output_dir) // "/" // trim(sim_name) // "_config.txt"
  call CFG_write(sim_cfg, trim(tmp_name))
  print *, "Settings written to ", trim(tmp_name)

  ! Initialize the transport coefficients
  call init_transport_coeff(sim_cfg)

  ! Set the initial conditions from the configuration
  applied_voltage = -domain_len * applied_fld
  call get_init_cond(sim_cfg, init_cond)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_fld
  mg%i_rhs        = i_rhs
  mg%i_eps        = i_eps

  ! The number of cycles at the lowest level
  mg%n_cycle_base = 8

  ! Routines to use for ...
  mg%sides_bc    => sides_bc_pot ! Filling ghost cell on physical boundaries
  mg%box_op      => mg2_auto_op
  mg%box_corr    => mg2_auto_corr
  mg%box_gsrb    => mg2_auto_gsrb

  ! This routine always needs to be called when using multigrid
  call mg2_init_mg(mg)

  output_cnt = 0          ! Number of output files written
  time       = 0          ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a2_loop_box(tree, set_init_cond)
     call compute_fld(tree, n_fmg_cycles, .true.)
     call a2_adjust_refinement(tree, set_ref_flags, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  if (photoi_enabled) &
       call set_photoionization(tree, photoi_eta, photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     dt      = get_max_dt(tree)

     if (dt < 1e-14) then
        print *, "dt getting too small, instability?"
        time = end_time + 1.0_dp
     end if

     ! Every dt_output, write output
     if (output_cnt * dt_output <= time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I6.6)") trim(sim_name) // "_", output_cnt
        fname_axis = trim(output_dir) // "/" // trim(fname) // "_axis.txt"
        fname_stats = trim(output_dir) // "/" // trim(sim_name) // ".txt"
     else
        write_out = .false.
     end if

     if (time > end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, n_steps_amr

        ! Copy previous solution
        call a2_tree_copy_cc(tree, i_elec, i_elec_old)
        call a2_tree_copy_cc(tree, i_pion, i_pion_old)

        ! Two forward Euler steps over dt
        do i = 1, 2
           time = time + dt

           ! First calculate fluxes
           call a2_loop_boxes_arg(tree, fluxes_koren, [dt], .true.)
           call a2_consistent_fluxes(tree, [f_elec])

           ! Update the solution
           call a2_loop_box_arg(tree, update_solution, [dt], .true.)

           ! Restrict the electron and ion densities to lower levels
           call a2_restrict_tree(tree, i_elec)
           call a2_restrict_tree(tree, i_pion)

           ! Fill ghost cells
           call a2_gc_sides(tree, i_elec, a2_sides_interp_lim, a2_bc_neumann)
           call a2_gc_sides(tree, i_pion, a2_sides_interp_lim, a2_bc_neumann)

           ! Compute new field on first iteration
           if (i == 1) call compute_fld(tree, n_fmg_cycles, .false.)
        end do

        time = time - dt        ! Go back one time step

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a2_loop_box(tree, average_dens)

        ! Compute field with new density
        call compute_fld(tree, n_fmg_cycles, .false.)
     end do

     if (write_out) then
        call a2_write_silo(tree, fname, &
             cc_names, output_cnt, time, dir=output_dir, &
             fc_names=["fld_r", "fld_z"], ixs_fc=[f_fld])
        call write_streamer_properties(tree, fname_stats, &
             fname_axis, output_cnt==1)
     end if

     call a2_adjust_refinement(tree, set_ref_flags, ref_info)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_fld(tree, n_fmg_cycles, .false.)

        ! This will every now-and-then clean up the data in the tree
        call a2_tidy_up(tree, 0.9_dp, 0.25_dp, 5000, .false.)
     end if

     if (photoi_enabled) &
          call set_photoionization(tree, photoi_eta, photoi_num_photons)
  end do

  call a2_destroy(tree)

contains

  ! Initialize the AMR tree
  subroutine init_tree(tree)
    type(a2_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: nb_list(4, 1) ! Neighbors of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = domain_len / box_size

    ! Initialize tree
    call a2_init(tree, box_size, n_var_cell, n_var_face, dr, &
         coarsen_to=4, n_boxes=n_boxes_init, coord=a5_cyl)

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...
    nb_list(:, id) = -1         ! And neighbors -1 (physical boundary)

    ! Create the base mesh
    call a2_set_base(tree, ix_list, nb_list)

  end subroutine init_tree

  ! Refinement function
  subroutine set_ref_flags(boxes, id, ref_flags)
    use m_geom
    type(box2_t), intent(in) :: boxes(:)
    integer, intent(in)      :: id
    integer, intent(inout)   :: ref_flags(:)
    integer                  :: nc, n
    real(dp)                 :: crv_phi, dr2, max_fld, max_dns
    real(dp)                 :: boxlen, dist, alpha, adx

    nc        = boxes(id)%n_cell
    dr2       = boxes(id)%dr**2
    crv_phi   = dr2 * maxval(abs(boxes(id)%cc(1:nc, 1:nc, i_rhs)))
    max_fld   = maxval(boxes(id)%cc(1:nc, 1:nc, i_fld))
    max_dns   = maxval(boxes(id)%cc(1:nc, 1:nc, i_elec))
    alpha     = LT_get_col(td_tbl, i_alpha, max_fld)
    adx       = boxes(id)%dr * alpha

    if (adx < 0.1_dp .and. boxes(id)%dr < 2.5e-5_dp) &
         ref_flags(id) = a5_rm_ref

    if (time < 5.0e-9_dp) then
       boxlen = boxes(id)%n_cell * boxes(id)%dr

       do n = 1, init_cond%n_cond
          dist = GM_dist_line(a2_r_center(boxes(id)), &
               init_cond%seed_r0(:, n), &
               init_cond%seed_r1(:, n), 2)
          if (dist - init_cond%seed_width(n) < boxlen &
               .and. boxes(id)%dr > 0.2_dp * init_cond%seed_width(n)) then
             ref_flags(id) = a5_do_ref
          end if
       end do
    end if

    if (adx > 0.75_dp .and. crv_phi > 0.1_dp) ref_flags(id) = a5_do_ref
  end subroutine set_ref_flags

  subroutine set_init_cond(box)
    use m_geom
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, n, nc
    real(dp)                    :: xy(2)
    real(dp)                    :: dens

    nc = box%n_cell
    box%cc(:, :, i_elec) = init_cond%bg_dens

    do j = 0, nc+1
       do i = 0, nc+1
          xy   = a2_r_cc(box, [i,j])

          do n = 1, init_cond%n_cond
             dens = init_cond%seed_dens(n) * &
                  GM_dens_line(xy, init_cond%seed_r0(:, n), &
                  init_cond%seed_r1(:, n), 2, &
                  init_cond%seed_width(n), &
                  init_cond%seed_falloff(n))
             box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + dens
          end do
       end do
    end do

    box%cc(:, :, i_pion) = box%cc(:, :, i_elec)
    box%cc(:, :, i_phi) = 0     ! Inital potential set to zero

    call set_box_eps(box)
  end subroutine set_init_cond

  subroutine set_box_eps(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell

    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j]) / domain_len

          if (xy(1) < 0.25_dp) then
             box%cc(i, j, i_eps) = epsilon_diel
          else
             box%cc(i, j, i_eps) = 1.0_dp
          end if
       end do
    end do
  end subroutine set_box_eps

  ! Get maximum time step based on e.g. CFL criteria
  real(dp) function get_max_dt(tree)
    type(a2_t), intent(in) :: tree
    real(dp), parameter    :: UC_eps0        = 8.8541878176d-12
    real(dp), parameter    :: UC_elem_charge = 1.6022d-19
    real(dp)               :: max_fld, min_fld, max_dns, dr_min
    real(dp)               :: mobility, diff_coeff, alpha, max_mobility
    real(dp)               :: dt_cfl, dt_dif, dt_drt, dt_alpha

    call a2_tree_max_cc(tree, i_fld, max_fld)
    call a2_tree_min_cc(tree, i_fld, min_fld)
    call a2_tree_max_cc(tree, i_elec, max_dns)

    dr_min       = a2_min_dr(tree)
    mobility     = LT_get_col(td_tbl, i_mobility, max_fld)
    max_mobility = LT_get_col(td_tbl, i_mobility, min_fld)
    diff_coeff   = LT_get_col(td_tbl, i_diffusion, max_fld)
    alpha        = LT_get_col(td_tbl, i_alpha, max_fld)

    ! CFL condition. Note there should be a factor sqrt(0.5), but instead we
    ! rely on the CFL number 0.5 at the bottom
    dt_cfl = dr_min / (mobility * max_fld)

    ! Diffusion condition
    dt_dif = 0.25_dp * dr_min**2 / diff_coeff

    ! Dielectric relaxation time
    dt_drt = UC_eps0 / (UC_elem_charge * max_mobility * &
         max(epsilon(1.0_dp), max_dns))

    ! Ionization limit
    dt_alpha =  1 / max(mobility * max_fld * alpha, epsilon(1.0_dp))

    get_max_dt = 0.5_dp * min(1/(1/dt_cfl + 1/dt_dif), dt_alpha, dt_max)
  end function get_max_dt

  ! Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_fld(tree, n_cycles, no_guess)
    use m_units_constants
    type(a2_t), intent(inout) :: tree
    integer, intent(in)       :: n_cycles
    logical, intent(in)       :: no_guess
    real(dp), parameter       :: fac = UC_elem_charge / UC_eps0
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
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

    applied_voltage = -domain_len * get_fld(time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg2_fas_fmg(tree, mg, .true., no_guess .and. i == 1)
    end do

    ! Compute field from potential
    call a2_loop_box(tree, fld_from_pot)

    ! Set the field norm also in ghost cells
    call a2_gc_sides(tree, i_fld, a2_sides_interp, a2_bc_neumann)
  end subroutine compute_fld

  real(dp) function get_fld(time)
    use m_units_constants
    real(dp), intent(in) :: time

    if (time > fld_mod_t0) then
       get_fld = applied_fld + (time - fld_mod_t0) * fld_lin_deriv + &
            fld_sin_amplitude * &
            sin((time - fld_mod_t0) * 2 * UC_pi * fld_sin_freq)
    else
       get_fld = applied_fld
    end if
  end function get_fld

  ! Compute electric field from electrical potential
  subroutine fld_from_pot(box)
    type(box2_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

    box%fx(:, :, f_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fy(:, :, f_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    ! Compute fields at the boundaries of the box, where eps can change (have to
    ! be careful that there is enough refinement)
    box%fx(1, :, f_fld) = 2 * inv_dr * &
         (box%cc(0, 1:nc, i_phi) - box%cc(1, 1:nc, i_phi)) * &
         box%cc(0, 1:nc, i_eps) / &
         (box%cc(1, 1:nc, i_eps) + box%cc(0, 1:nc, i_eps))
    box%fx(nc+1, :, f_fld) = 2 * inv_dr * &
         (box%cc(nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, i_phi)) * &
         box%cc(nc+1, 1:nc, i_eps) / &
         (box%cc(nc+1, 1:nc, i_eps) + box%cc(nc, 1:nc, i_eps))
    box%fy(:, 1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, 0, i_phi) - box%cc(1:nc, 1, i_phi)) * &
         box%cc(1:nc, 0, i_eps) / &
         (box%cc(1:nc, 1, i_eps) + box%cc(1:nc, 0, i_eps))
    box%fy(:, nc+1, f_fld) = 2 * inv_dr * &
         (box%cc(1:nc, nc, i_phi) - box%cc(1:nc, nc+1, i_phi)) * &
         box%cc(1:nc, nc+1, i_eps) / &
         (box%cc(1:nc, nc+1, i_eps) + box%cc(1:nc, nc, i_eps))

    box%cc(1:nc, 1:nc, i_fld) = sqrt(&
         0.25_dp * (box%fx(1:nc, 1:nc, f_fld) + box%fx(2:nc+1, 1:nc, f_fld))**2 + &
         0.25_dp * (box%fy(1:nc, 1:nc, f_fld) + box%fy(1:nc, 2:nc+1, f_fld))**2)
  end subroutine fld_from_pot

  ! This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_pot(boxes, id, nb, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv
    integer                     :: nc

    nc = boxes(id)%n_cell

    select case (nb)
    case (a2_nb_lx)             ! Neumann
       boxes(id)%cc(0, 1:nc, iv) = boxes(id)%cc(1, 1:nc, iv)
    case (a2_nb_hx)             ! Neumann
       boxes(id)%cc(nc+1, 1:nc, iv) = boxes(id)%cc(nc, 1:nc, iv)
    case (a2_nb_ly)             ! Grounded
       boxes(id)%cc(1:nc, 0, iv) = -boxes(id)%cc(1:nc, 1, iv)
    case (a2_nb_hy)             ! Applied voltage
       boxes(id)%cc(:, nc+1, iv) = 2 * applied_voltage &
            - boxes(id)%cc(:, nc, iv)
    end select
  end subroutine sides_bc_pot

  !> Modified implementation of Koren limiter, to avoid division and the min/max
  !> functions, which can be problematic / expensive. In most literature, you
  !> have r = a / b (ratio of gradients). Then the limiter phi(r) is multiplied
  !> with b. With this implementation, you get phi(r) * b
  elemental function koren_mlim(a, b) result(bphi)
    real(dp), intent(in) :: a  !< Density gradient (numerator)
    real(dp), intent(in) :: b  !< Density gradient (denominator)
    real(dp), parameter  :: sixth = 1/6.0_dp
    real(dp)             :: bphi, aa, ab

    aa = a * a
    ab = a * b

    if (ab <= 0) then
       ! a and b have different sign or one of them is zero, so r is either 0,
       ! inf or negative (special case a == b == 0 is ignored)
       bphi = 0
    else if (aa <= 0.25_dp * ab) then
       ! 0 < a/b <= 1/4, limiter has value a/b
       bphi = a
    else if (aa <= 2.5_dp * ab) then
       ! 1/4 < a/b <= 2.5, limiter has value (1+2*a/b)/6
       bphi = sixth * (b + 2*a)
    else
       ! (1+2*a/b)/6 >= 1, limiter has value 1
       bphi = b
    end if
  end function koren_mlim

  ! Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id, dt_vec)
    use m_units_constants
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: dt_vec(:)
    real(dp)                    :: fac, inv_dr, tmp, gradp, gradc, gradn
    real(dp)                    :: mobility, diff_coeff, v_drift
    real(dp)                    :: fld, fld_avg
    real(dp)                    :: gc_data(boxes(id)%n_cell, a2_num_neighbors)
    integer                     :: i, j, nc
    type(LT_loc_t) :: loc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr
    fac    = -0.8_dp * UC_eps0 / (UC_elem_charge * dt_vec(1))

    call a2_gc2_box_sides(boxes, id, i_elec, a2_sides2_prolong1, &
         sides_bc2_dens, gc_data, nc)

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          fld_avg   = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i-1, j, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fx(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i-1, j, i_elec)

          if (v_drift < 0.0_dp) then
             if (i == nc+1) then
                tmp = gc_data(j, a2_nb_hx)
             else
                tmp = boxes(id)%cc(i+1, j, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
             if (boxes(id)%fx(i, j, f_elec) < fac * fld) &
                  boxes(id)%fx(i, j, f_elec) = fac * fld
          else                  ! v_drift > 0
             if (i == 1) then
                tmp = gc_data(j, a2_nb_lx)
             else
                tmp = boxes(id)%cc(i-2, j, i_elec)
             end if
             gradp = boxes(id)%cc(i-1, j, i_elec) - tmp
             boxes(id)%fx(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i-1, j, i_elec) + koren_mlim(gradc, gradp))
             if (boxes(id)%fx(i, j, f_elec) > fac * fld) &
                  boxes(id)%fx(i, j, f_elec) = fac * fld
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fx(i, j, f_elec) = boxes(id)%fx(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          fld_avg    = 0.5_dp * (boxes(id)%cc(i, j, i_fld) + &
               boxes(id)%cc(i, j-1, i_fld))
          loc        = LT_get_loc(td_tbl, fld_avg)
          mobility   = LT_get_col_at_loc(td_tbl, i_mobility, loc)
          diff_coeff = LT_get_col_at_loc(td_tbl, i_diffusion, loc)
          fld        = boxes(id)%fy(i, j, f_fld)
          v_drift    = -mobility * fld
          gradc      = boxes(id)%cc(i, j, i_elec) - boxes(id)%cc(i, j-1, i_elec)

          if (v_drift < 0.0_dp) then
             if (j == nc+1) then
                tmp = gc_data(i, a2_nb_hy)
             else
                tmp = boxes(id)%cc(i, j+1, i_elec)
             end if
             gradn = tmp - boxes(id)%cc(i, j, i_elec)
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j, i_elec) - koren_mlim(gradc, gradn))
             if (boxes(id)%fy(i, j, f_elec) < fac * fld) &
                  boxes(id)%fy(i, j, f_elec) = fac * fld
          else                  ! v_drift > 0
             if (j == 1) then
                tmp = gc_data(i, a2_nb_ly)
             else
                tmp = boxes(id)%cc(i, j-2, i_elec)
             end if
             gradp = boxes(id)%cc(i, j-1, i_elec) - tmp
             boxes(id)%fy(i, j, f_elec) = v_drift * &
                  (boxes(id)%cc(i, j-1, i_elec) + koren_mlim(gradc, gradp))
             if (boxes(id)%fy(i, j, f_elec) > fac * fld) &
                  boxes(id)%fy(i, j, f_elec) = fac * fld
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be
          ! scaled by 1/dx
          boxes(id)%fy(i, j, f_elec) = boxes(id)%fy(i, j, f_elec) - &
               diff_coeff * gradc * inv_dr
       end do
    end do

  end subroutine fluxes_koren

  ! Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_dens(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_elec) = 0.5_dp * (box%cc(:, :, i_elec) + box%cc(:, :, i_elec_old))
    box%cc(:, :, i_pion) = 0.5_dp * (box%cc(:, :, i_pion) + box%cc(:, :, i_pion_old))
  end subroutine average_dens

  ! Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr, src, sflux, fld
    real(dp)                    :: alpha, eta, dflux(2), rfac(2)
    integer                     :: i, j, nc, ioff
    type(LT_loc_t)              :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
    ioff   = (box%ix(1)-1) * nc

    do j = 1, nc
       do i = 1, nc
          ! Weighting of flux contribution for cylindrical coordinates
          rfac = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)

          fld      = box%cc(i,j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          eta      = LT_get_col_at_loc(td_tbl, i_eta, loc)

          ! Set source term equal to ||flux|| * (alpha - eta)
          dflux(1) = box%fx(i, j, f_elec) + box%fx(i+1, j, f_elec)
          dflux(2) = box%fy(i, j, f_elec) + box%fy(i, j+1, f_elec)
          src = 0.5_dp * norm2(dflux) * (alpha - eta)

          if (photoi_enabled) &
               src = src + box%cc(i,j, i_pho)

          ! Contribution of flux
          sflux = (box%fy(i, j, f_elec) - box%fy(i, j+1, f_elec) + &
               rfac(1) * box%fx(i, j, f_elec) - &
               rfac(2) * box%fx(i+1, j, f_elec)) * inv_dr

          box%cc(i, j, i_elec) = box%cc(i, j, i_elec) + (src + sflux) * dt(1)
          box%cc(i, j, i_pion) = box%cc(i, j, i_pion) + src * dt(1)
       end do
    end do
  end subroutine update_solution

  subroutine set_photoionization(tree, eta, num_photons)
    use m_units_constants

    type(a2_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0D0 * UC_torr_to_bar
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a2_loop_box_arg(tree, set_photoi_rate, [eta * quench_fac], .true.)

    call PH_set_src_2d(tree, photoi_tbl, sim_rng, num_photons, &
         i_pho, i_pho, 0.6_dp, .false., .true., 0.05e-3_dp)

  end subroutine set_photoionization

  subroutine set_photoi_rate(box, coeff)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: i, j, nc
    real(dp)                    :: fld, alpha, mobility, dr, tmp
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do j = 1, nc
       do i = 1, nc
          dr       = box%dr
          fld      = box%cc(i, j, i_fld)
          loc      = LT_get_loc(td_tbl, fld)
          alpha    = LT_get_col_at_loc(td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(td_tbl, i_mobility, loc)

          tmp = fld * mobility * alpha * box%cc(i, j, i_elec) * coeff(1)
          if (tmp < 0) tmp = 0
          box%cc(i, j, i_pho) = tmp
       end do
    end do
  end subroutine set_photoi_rate

  ! For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id

    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%max_lvl
       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_prolong1_to(tree%boxes, id, i_elec)
          call a2_prolong1_to(tree%boxes, id, i_pion)
          call a2_prolong1_to(tree%boxes, id, i_phi)
          call set_box_eps(tree%boxes(id))
       end do
       !$omp end do

       !$omp do
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a2_gc_box_sides(tree%boxes, id, i_elec, &
               a2_sides_interp_lim, a2_bc_neumann)
          call a2_gc_box_sides(tree%boxes, id, i_pion, &
               a2_sides_interp_lim, a2_bc_neumann)
          call a2_gc_box_sides(tree%boxes, id, i_phi, &
               mg2_sides_rb, sides_bc_pot)
       end do
       !$omp end do
    end do
    !$omp end parallel
  end subroutine prolong_to_new_boxes

  ! This fills a second layer of ghost cells near physical boundaries for the
  ! electron density
  subroutine sides_bc2_dens(boxes, id, nb, iv, bc_side, nc)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, nb, iv, nc
    real(dp), intent(out)       :: bc_side(nc)

    select case (nb)
    case (a2_nb_lx)
       ! Neumann zero
       bc_side = boxes(id)%cc(2, 1:nc, iv)
    case (a2_nb_hx)
       ! Neumann zero
       bc_side = boxes(id)%cc(nc-1, 1:nc, iv)
    case (a2_nb_ly)
       ! Neumann zero
       bc_side = boxes(id)%cc(1:nc, 2, iv)
    case (a2_nb_hy)
       ! Neumann zero
       bc_side = boxes(id)%cc(1:nc, nc-1, iv)
    end select
  end subroutine sides_bc2_dens

  subroutine get_init_cond(cfg, cond)
    type(CFG_t), intent(in)        :: cfg
    type(initcnd_t), intent(inout) :: cond
    integer                        :: n_cond, varsize
    real(dp)                       :: dlen
    real(dp), allocatable          :: tmp_vec(:)

    call CFG_get(cfg, "bg_dens", cond%bg_dens)
    call CFG_get(cfg, "domain_len", dlen)

    call CFG_get_size(cfg, "seed_dens", n_cond)

    call CFG_get_size(cfg, "seed_rel_r0", varsize)
    if (varsize /= 2 * n_cond) &
         stop "seed_... variables have incompatible size"

    call CFG_get_size(cfg, "seed_rel_r1", varsize)
    if (varsize /= 2 * n_cond) &
         stop "seed_... variables have incompatible size"

    call CFG_get_size(cfg, "seed_width", varsize)
    if (varsize /= n_cond) &
         stop "seed_... variables have incompatible size"

    cond%n_cond = n_cond
    allocate(cond%seed_dens(n_cond))
    allocate(cond%seed_r0(2, n_cond))
    allocate(cond%seed_r1(2, n_cond))
    allocate(cond%seed_width(n_cond))
    allocate(cond%seed_falloff(n_cond))

    allocate(tmp_vec(2 * n_cond))
    call CFG_get(cfg, "seed_rel_r0", tmp_vec)
    cond%seed_r0 = dlen * reshape(tmp_vec, [2, n_cond])
    call CFG_get(cfg, "seed_rel_r1", tmp_vec)
    cond%seed_r1 = dlen * reshape(tmp_vec, [2, n_cond])

    call CFG_get(cfg, "seed_dens", cond%seed_dens)
    call CFG_get(cfg, "seed_width", cond%seed_width)
    call CFG_get(cfg, "seed_falloff", cond%seed_falloff)
  end subroutine get_init_cond

  subroutine create_cfg(cfg)
    type(CFG_t), intent(inout) :: cfg

    call CFG_add(cfg, "end_time", 10.0d-9, &
         "The desired endtime in seconds of the simulation")
    call CFG_add(cfg, "sim_name", "sim", &
         "The name of the simulation")
    call CFG_add(cfg, "output_dir", "", &
         "Directory where the output should be written")
    call CFG_add(cfg, "box_size", 8, &
         "The number of grid cells per coordinate in a box")
    call CFG_add(cfg, "domain_len", 32e-3_dp, &
         "The length of the (square) domain")
    call CFG_add(cfg, "gas_name", "N2", &
         "The name of the gas mixture used")
    call CFG_add(cfg, "gas_pressure", 1.0_dp, &
         "The gas pressure in bar (used for photoionization)")
    call CFG_add(cfg, "applied_fld", 1.0d7, &
         "The applied electric field")
    call CFG_add(cfg, "epsilon_diel", 1.5_dp, &
         "The dielectric constant of the dielectric")

    call CFG_add(cfg, "fld_mod_t0", 1.0e99_dp, &
         "Modify electric field after this time")
    call CFG_add(cfg, "fld_sin_amplitude", 0.0_dp, &
         "Amplitude of sinusoidal modification")
    call CFG_add(cfg, "fld_sin_freq", 0.2e9_dp, &
         "Frequency of sinusoidal modification")
    call CFG_add(cfg, "fld_lin_deriv", 0.0_dp, &
         "Linear derivative of field")

    call CFG_add(cfg, "bg_dens", 1.0d12, &
         "The background ion and electron density in 1/m^3")
    call CFG_add(cfg, "seed_dens", [5.0d19], &
         "Initial density of the seed", .true.)
    call CFG_add(cfg, "seed_rel_r0", [0.5d0, 0.4d0], &
         "The relative start position of the initial seed", .true.)
    call CFG_add(cfg, "seed_rel_r1", [0.5d0, 0.6d0], &
         "The relative end position of the initial seed", .true.)
    call CFG_add(cfg, "seed_width", [0.5d-3], &
         "Seed width", .true.)
    call CFG_add(cfg, "seed_falloff", [1], &
         "Fallof type for seed, see m_geom.f90", .true.)

    call CFG_add(cfg, "dt_output", 1.0d-10, &
         "The timestep for writing output")
    call CFG_add(cfg, "dt_max", 1.0d-11, &
         "The maximum timestep")
    call CFG_add(cfg, "num_steps_amr", 2, &
         "The number of steps after which the mesh is updated")

    call CFG_add(cfg, "photoi_enabled", .true., &
         "Whether photoionization is enabled")
    call CFG_add(cfg, "photoi_frac_O2", 0.2_dp, &
         "Fraction of oxygen")
    call CFG_add(cfg, "photoi_eta", 0.05_dp, &
         "Photoionization efficiency factor")
    call CFG_add(cfg, "photoi_num_photons", 50*1000, &
         "Number of discrete photons to use for photoionization")

    call CFG_add(cfg, "input_file", "transport_data_file.txt", &
         "Input file with transport data")
    call CFG_add(cfg, "lkptbl_size", 1000, &
         "The transport data table size in the fluid model")
    call CFG_add(cfg, "lkptbl_max_fld", 3.0d7, &
         "The maximum electric field in the fluid model coefficients")
    call CFG_add(cfg, "td_mobility_name", "efield[V/m]_vs_mu[m2/Vs]", &
         "The name of the mobility coefficient")
    call CFG_add(cfg, "td_diffusion_name", "efield[V/m]_vs_dif[m2/s]", &
         "The name of the diffusion coefficient")
    call CFG_add(cfg, "td_alpha_name", "efield[V/m]_vs_alpha[1/m]", &
         "The name of the eff. ionization coeff.")
    call CFG_add(cfg, "td_eta_name", "efield[V/m]_vs_eta[1/m]", &
         "The name of the eff. attachment coeff.")

    call CFG_add(cfg, "td_alpha_fac", 1.0_dp, &
         "Modify alpha by this factor")
    call CFG_add(cfg, "td_eta_fac", 1.0_dp, &
         "Modify eta by this factor")
    call CFG_add(cfg, "td_mobility_fac", 1.0_dp, &
         "Modify mobility by this factor")
    call CFG_add(cfg, "td_diffusion_fac", 1.0_dp, &
         "Modify diffusion by this factor")
  end subroutine create_cfg

  subroutine init_transport_coeff(cfg)
    use m_transport_data
    use m_config

    type(CFG_t), intent(in) :: cfg
    character(len=name_len) :: input_file, gas_name
    integer                 :: table_size
    real(dp)                :: max_fld, alpha_fac, eta_fac
    real(dp)                :: mobility_fac, diffusion_fac
    real(dp), allocatable   :: x_data(:), y_data(:)
    character(len=name_len) :: data_name

    call CFG_get(cfg, "input_file", input_file)
    call CFG_get(cfg, "gas_name", gas_name)

    call CFG_get(cfg, "lkptbl_size", table_size)
    call CFG_get(cfg, "lkptbl_max_fld", max_fld)

    call CFG_get(cfg, "td_alpha_fac", alpha_fac)
    call CFG_get(cfg, "td_eta_fac", eta_fac)
    call CFG_get(cfg, "td_mobility_fac", mobility_fac)
    call CFG_get(cfg, "td_diffusion_fac", diffusion_fac)

    ! Create a lookup table for the model coefficients
    td_tbl = LT_create(0.0_dp, max_fld, table_size, n_var_td)

    ! Fill table with data
    call CFG_get(cfg, "td_mobility_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * mobility_fac
    call LT_set_col(td_tbl, i_mobility, x_data, y_data)

    call CFG_get(cfg, "td_diffusion_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * diffusion_fac
    call LT_set_col(td_tbl, i_diffusion, x_data, y_data)

    call CFG_get(cfg, "td_alpha_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * alpha_fac
    call LT_set_col(td_tbl, i_alpha, x_data, y_data)

    call CFG_get(cfg, "td_eta_name", data_name)
    call TD_get_td_from_file(input_file, gas_name, &
         trim(data_name), x_data, y_data)
    y_data = y_data * eta_fac
    call LT_set_col(td_tbl, i_eta, x_data, y_data)

    ! Create table for photoionization
    call CFG_get(cfg, "gas_pressure", gas_pressure)
    call CFG_get(cfg, "photoi_enabled", photoi_enabled)
    call CFG_get(cfg, "photoi_frac_O2", photoi_frac_O2)
    call CFG_get(cfg, "photoi_eta", photoi_eta)
    call CFG_get(cfg, "photoi_num_photons", photoi_num_photons)

    if (photoi_enabled) then
       call PH_get_tbl_air(photoi_tbl, photoi_frac_O2 * gas_pressure, &
            2 * domain_len)
    end if

  end subroutine init_transport_coeff

  subroutine write_streamer_properties(tree, fname_stats, fname_axis, first_time)
    type(a2_t), intent(in)       :: tree
    character(len=*), intent(in) :: fname_axis, fname_stats
    logical, intent(in)          :: first_time

    real(dp)                     :: fld_z, fld_r, radius, radius_z, height
    real(dp)                     :: edens, phi, bg_fld
    real(dp), allocatable        :: axis_data(:,:)
    integer                      :: n
    integer, parameter           :: unit_1 = 777, unit_2 = 778

    call get_streamer_properties(tree, height, fld_z, radius, radius_z, &
         fld_r, edens, phi, bg_fld)
    call get_cc_axis(tree, [i_elec, i_pion], [f_fld, f_elec], axis_data)

    if (first_time) then
       open(unit_1, file=trim(fname_stats), action="write")
       write(unit_1, *) "#time, height, fld_z, radius, radius_z, ", &
            "fld_r, edens, phi, bg_fld"
       close(unit_1)
    else
       open(unit_1, file=trim(fname_stats), action="write", &
            position="append")
       write(unit_1, *) time, height, fld_z, radius, radius_z, &
            fld_r, edens, phi, bg_fld
       close(unit_1)
    end if

    open(unit_2, file=trim(fname_axis), action="write")
    write(unit_1, *) "#z, n_e, n_i, fld_z, flux_z"
    do n = 1, size(axis_data, 2)
       write(unit_2, *) axis_data(:, n)
    end do
    close(unit_2)
    deallocate(axis_data)

    print *, "Written ", trim(fname_axis)
    if (.not. first_time .and. height > 0.9_dp * domain_len) &
         stop "Simulation has reached boundary"
  end subroutine write_streamer_properties

  subroutine get_streamer_properties(tree, height, fld_z, &
       radius, radius_z, fld_r, edens, phi, bg_fld)
    type(a2_t), intent(in) :: tree
    real(dp), intent(out) :: fld_z, fld_r, radius, radius_z, height
    real(dp), intent(out) :: edens, phi, bg_fld

    real(dp)       :: rz(2)
    type(a2_loc_t) :: loc_ez, loc_er, loc_dens
    integer        :: id, ix(2)

    call a2_reduction_loc(tree, box_fld_z, reduce_max, &
         -1.0e99_dp, fld_z, loc_ez)
    call a2_reduction_loc(tree, box_fld_r, reduce_max, &
         -1.0e99_dp, fld_r, loc_er)

    ! Radius of streamer is defined as location of maximum r-field
    rz       = a2_r_cc(tree%boxes(loc_er%id), loc_er%ix)
    radius   = rz(1)
    radius_z = rz(2)

    ! Get electron density and potential at location of radius
    loc_dens = a2_get_loc(tree, [0.0_dp, rz(2)])
    id       = loc_dens%id
    ix       = loc_dens%ix
    edens    = tree%boxes(id)%cc(ix(1), ix(2), i_elec)
    ! Set phi to potential difference
    phi      = tree%boxes(id)%cc(ix(1), ix(2), i_phi)
    phi      = phi - (rz(2)/domain_len) * applied_voltage

    ! Height of streamer
    rz     = a2_r_cc(tree%boxes(loc_ez%id), loc_ez%ix)
    height = rz(2)

    bg_fld = get_fld(time)
  end subroutine get_streamer_properties

  real(dp) function reduce_max(a, b)
    real(dp), intent(in) :: a, b
    reduce_max = max(a,b)
  end function reduce_max

  subroutine box_fld_r(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: nc

    nc = box%n_cell
    ix = maxloc(abs(box%fx(1:nc, :, f_fld) + box%fx(2:nc+1, :, f_fld)))
    val = 0.5_dp * abs(box%fx(ix(1), ix(2), f_fld) + &
         box%fx(ix(1)+1, ix(2), f_fld))
  end subroutine box_fld_r

  subroutine box_fld_z(box, val, ix)
    type(box2_t), intent(in) :: box
    real(dp), intent(out)    :: val
    integer, intent(out)     :: ix(2)
    integer                  :: nc

    nc = box%n_cell
    ix = maxloc(abs(box%fy(:, 1:nc, f_fld) + box%fy(:, 2:nc+1, f_fld)))
    val = 0.5_dp * abs(box%fy(ix(1), ix(2), f_fld) + &
         box%fy(ix(1), ix(2)+1, f_fld))
  end subroutine box_fld_z

  subroutine get_cc_axis(tree, ixs_cc, ixs_fc, axis_data)
    type(a2_t), intent(in)               :: tree
    integer, intent(in)                  :: ixs_cc(:), ixs_fc(:)
    real(dp), allocatable, intent(inout) :: axis_data(:, :)

    type(a2_loc_t)                       :: loc
    real(dp)                             :: z, box_z, box_dz
    integer                              :: i, id, nc, cnt, n_cc, n_fc

    n_cc = size(ixs_cc)
    n_fc = size(ixs_fc)
    nc   = tree%n_cell
    z    = 0
    cnt  = 0

    ! Determine how many boxes lie on the axis
    do
       loc = a2_get_loc(tree, [0.0_dp, z])
       if (loc%id == -1) exit

       cnt    = cnt + nc
       id     = loc%id
       box_z  = tree%boxes(id)%r_min(2)
       box_dz = tree%boxes(id)%dr
       z      = box_z + (nc+1) * box_dz
    end do

    ! Now store the actual axis data
    allocate(axis_data(n_cc+n_fc+1, cnt))
    cnt = 0
    z   = 0

    do
       loc = a2_get_loc(tree, [0.0_dp, z])
       if (loc%id == -1) exit

       id     = loc%id
       box_z  = tree%boxes(id)%r_min(2)
       box_dz = tree%boxes(id)%dr

       axis_data(1, cnt+1:cnt+nc) = &
            box_z + [((i-0.5_dp) * box_dz, i = 1, nc)]
       axis_data(2:n_cc+1, cnt+1:cnt+nc) = &
            transpose(tree%boxes(id)%cc(1, 1:nc, ixs_cc))
       axis_data(n_cc+2:, cnt+1:cnt+nc) = transpose( &
            0.5_dp * (tree%boxes(id)%fy(1, 1:nc, ixs_fc) + &
            tree%boxes(id)%fy(1, 2:nc+1, ixs_fc)))

       cnt    = cnt + nc
       z      = box_z + (nc+1) * box_dz
    end do
  end subroutine get_cc_axis

end program streamer_cyl
