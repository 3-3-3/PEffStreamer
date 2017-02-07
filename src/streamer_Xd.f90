#include "../afivo/src/cpp_macros_$Dd.h"
!> Program to perform $Dd streamer simulations in Cartesian and cylindrical coordinates
program streamer_$Dd

  use m_a$D_types
  use m_a$D_core
  use m_a$D_ghostcell
  use m_a$D_utils
  use m_a$D_restrict
  use m_a$D_multigrid
  use m_a$D_output
  use m_write_silo
  use m_streamer

  implicit none

  integer                :: i, n
  character(len=ST_slen) :: fname
  logical                :: write_out

  type(a$D_t)             :: tree ! This contains the full grid information
  type(mg$D_t)            :: mg   ! Multigrid option struct
  type(ref_info_t)       :: ref_info

  call ST_create_config()
  call ST_read_config_files()
  call ST_load_config()

  ! Initialize the transport coefficients
  call ST_load_transport_data()

  ! Set the initial conditions from the configuration
  call ST_get_init_cond($D)

  ! Initialize the tree (which contains all the mesh information)
  call init_tree(tree)

  ! Set the multigrid options. First define the variables to use
  mg%i_phi        = i_phi
  mg%i_tmp        = i_electric_fld
  mg%i_rhs        = i_rhs

  ! Routines to use for ...
  mg%sides_bc    => sides_bc_potential ! Filling ghost cell on physical boundaries

  ! This routine always needs to be called when using multigrid
  call mg$D_init_mg(mg)

  ST_out_cnt = 0 ! Number of output files written
  ST_time    = 0 ! Simulation time (all times are in s)

  ! Set up the initial conditions
  do
     call a$D_loop_box(tree, set_initial_condition)
     call compute_electric_field(tree, n_fmg_cycles, .false.)
     call a$D_adjust_refinement(tree, refine_routine, ref_info, 4)
     if (ref_info%n_add == 0) exit
  end do

  call a$D_print_info(tree)

  if (ST_photoi_enabled) &
       call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons)

  do
     ! Get a new time step, which is at most dt_amr
     call a$D_reduction_vec(tree, get_max_dt, ST_get_min, &
          [ST_dt_max, ST_dt_max, ST_dt_max], ST_dt_vec, ST_dt_num_cond)
     ST_dt = minval(ST_dt_vec)

     if (ST_dt < 1e-14) then
        print *, "ST_dt getting too small, instability?"
        ST_time = ST_end_time + 1.0_dp
     end if

     ! Every ST_dt_out, write output
     if (ST_out_cnt * ST_dt_out <= ST_time) then
        write_out = .true.
        ST_out_cnt = ST_out_cnt + 1
        write(fname, "(A,I6.6)") trim(ST_simulation_name) // "_", ST_out_cnt
     else
        write_out = .false.
     end if

     if (ST_time > ST_end_time) exit

     ! We perform n_steps between mesh-refinements
     do n = 1, ST_refine_per_steps

        if (ST_photoi_enabled) &
             call set_photoionization(tree, ST_photoi_eta, ST_photoi_num_photons, ST_dt)

        ! Copy previous solution
        call a$D_tree_copy_cc(tree, i_electron, i_electron_old)
        call a$D_tree_copy_cc(tree, i_pos_ion, i_pos_ion_old)

        ! Two forward Euler steps over ST_dt
        do i = 1, 2
           ST_time = ST_time + ST_dt

           ! First calculate fluxes
           call a$D_loop_boxes(tree, fluxes_koren, .true.)
           call a$D_consistent_fluxes(tree, [flux_elec])

           ! Update the solution
           call a$D_loop_box_arg(tree, update_solution, [ST_dt], .true.)

           ! Compute new field on first iteration
           if (i == 1) call compute_electric_field(tree, n_fmg_cycles, .true.)
        end do

        ST_time = ST_time - ST_dt        ! Go back one time step

        ! Take average of phi_old and phi (explicit trapezoidal rule)
        call a$D_loop_box(tree, average_density)

        ! Compute field with new density
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end do

     ! Restrict the electron and ion densities before refinement
     call a$D_restrict_tree(tree, i_electron)
     call a$D_restrict_tree(tree, i_pos_ion)

     ! Fill ghost cells before refinement
     call a$D_gc_tree(tree, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero)
     call a$D_gc_tree(tree, i_pos_ion, a$D_gc_interp_lim, a$D_bc_neumann_zero)

     ! if (write_out) call a$D_write_silo(tree, fname, ST_out_cnt, &
          ! ST_time, dir=ST_output_dir)

     call a$D_adjust_refinement(tree, refine_routine, ref_info, 4)

     if (ref_info%n_add > 0 .or. ref_info%n_rm > 0) then
        ! For boxes which just have been refined, set data on their children
        call prolong_to_new_boxes(tree, ref_info)

        ! Compute the field on the new mesh
        call compute_electric_field(tree, n_fmg_cycles, .true.)
     end if

  end do

  call a$D_destroy(tree)

contains

  !> Initialize the AMR tree
  subroutine init_tree(tree)
    type(a$D_t), intent(inout) :: tree

    ! Variables used below to initialize tree
    real(dp)                  :: dr
    integer                   :: id
    integer                   :: ix_list(2, 1) ! Spatial indices of initial boxes
    integer                   :: n_boxes_init = 1000

    dr = ST_domain_len / ST_box_size

    ! Initialize tree
    if (ST_cylindrical) then
       call a$D_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, coord=af_cyl, &
            cc_names=ST_cc_names)
    else
       call a$D_init(tree, ST_box_size, n_var_cell, n_var_face, dr, &
            coarsen_to=2, n_boxes=n_boxes_init, cc_names=ST_cc_names)
    end if

    ! Set up geometry
    id             = 1          ! One box ...
    ix_list(:, id) = [1,1]      ! With index 1,1 ...

    ! Create the base mesh
    call a$D_set_base(tree, 1, ix_list)

  end subroutine init_tree

  ! This routine sets the cell refinement flags for box
  subroutine refine_routine(box, cell_flags)
    use m_geometry
    type(box$D_t), intent(in) :: box
    ! Refinement flags for the cells of the box
    integer, intent(out)     :: &
         cell_flags(DTIMES(box%n_cell))
    integer                  :: IJK, n, nc
    real(dp)                 :: cphi, dx, dx2
    real(dp)                 :: alpha, adx, fld
    real(dp)                 :: dist

    nc      = box%n_cell
    dx      = box%dr
    dx2     = box%dr**2

    do KJI_DO(1,nc)
          fld   = box%cc(IJK, i_electric_fld)
          alpha = LT_get_col(ST_td_tbl, i_alpha, fld)
          ! The refinement is based on the ionization length
          adx   = box%dr * alpha

          ! The refinement is also based on the intensity of the source term.
          ! Here we estimate the curvature of phi (given by dx**2 *
          ! Laplacian(phi))
          cphi = dx2 * abs(box%cc(IJK, i_rhs))

          if (adx / ST_refine_adx + cphi / ST_refine_cphi > 1) then
             cell_flags(IJK) = af_do_ref
          else if (adx < 0.125_dp * ST_refine_adx .and. &
               cphi < 0.0625_dp * ST_refine_cphi &
               .and. dx < ST_derefine_dx) then
             cell_flags(IJK) = af_rm_ref
          else
             cell_flags(IJK) = af_keep_ref
          end if

          ! Refine around the initial conditions
          if (ST_time < ST_refine_init_time) then
             do n = 1, ST_init_cond%n_cond
                dist = GM_dist_line(a$D_r_cc(box, [IJK]), &
                     ST_init_cond%seed_r0(:, n), &
                     ST_init_cond%seed_r1(:, n), 2)
                if (dist - ST_init_cond%seed_width(n) < 2 * dx &
                     .and. box%dr > ST_refine_init_fac * &
                     ST_init_cond%seed_width(n)) then
                   cell_flags(IJK) = af_do_ref
                end if
             end do
          end if

    end do; CLOSE_DO

    ! Make sure we don't have or get a too fine or too coarse grid
    if (dx > ST_refine_max_dx) then
       cell_flags = af_do_ref
    else if (dx < 2 * ST_refine_min_dx) then
       where(cell_flags == af_do_ref) cell_flags = af_keep_ref
    end if

  end subroutine refine_routine

  !> Sets the initial condition
  subroutine set_initial_condition(box)
    use m_geometry
    type(box$D_t), intent(inout) :: box
    integer                     :: IJK, n, nc
    real(dp)                    :: rr($D)
    real(dp)                    :: density

    nc = box%n_cell
    box%cc(DTIMES(:), i_electron) = ST_init_cond%background_density
    box%cc(DTIMES(:), i_pos_ion)  = ST_init_cond%background_density
    box%cc(DTIMES(:), i_phi)      = 0 ! Inital potential set to zero

    do KJI_DO(0,nc+1)
          rr   = a$D_r_cc(box, [IJK])

          do n = 1, ST_init_cond%n_cond
             density = ST_init_cond%seed_density(n) * &
                  GM_density_line(rr, ST_init_cond%seed_r0(:, n), &
                  ST_init_cond%seed_r1(:, n), $D, &
                  ST_init_cond%seed_width(n), &
                  ST_init_cond%seed_falloff(n))

             ! Add electrons and/or ions depending on the seed charge type
             ! (positive, negative or neutral)
             if (ST_init_cond%seed_charge_type(n) <= 0) then
                box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + density
             end if

             if (ST_init_cond%seed_charge_type(n) >= 0) then
                box%cc(IJK, i_pos_ion) = box%cc(IJK, i_pos_ion) + density
             end if
          end do
    end do; CLOSE_DO

  end subroutine set_initial_condition

  !> Get maximum time step based on e.g. CFL criteria
  function get_max_dt(box, n_cond) result(dt_vec)
    use m_units_constants
    type(box$D_t), intent(in) :: box
    integer, intent(in)      :: n_cond
    integer                  :: IJK, nc
    real(dp)                 :: fld($D), fld_norm, mobility, diffusion_c
    real(dp)                 :: dt_vec(n_cond)

    nc = box%n_cell
    dt_vec = ST_dt_max

    do KJI_DO(1,nc)
#if $D == 2
       fld(1) = 0.5_dp * (box%fc(IJK, 1, electric_fld) + &
            box%fc(i+1, j, 1, electric_fld))
       fld(2) = 0.5_dp * (box%fc(IJK, 2, electric_fld) + &
            box%fc(i, j+1, 2, electric_fld))
#elif $D == 3
       fld(1) = 0.5_dp * (box%fc(IJK, 1, electric_fld) + &
            box%fc(i+1, j, k, 1, electric_fld))
       fld(2) = 0.5_dp * (box%fc(IJK, 2, electric_fld) + &
            box%fc(i, j+1, k, 2, electric_fld))
       fld(3) = 0.5_dp * (box%fc(IJK, 3, electric_fld) + &
            box%fc(i, j, k+1, 3, electric_fld))
#endif


       fld_norm = box%cc(IJK, i_electric_fld)
       mobility = LT_get_col(ST_td_tbl, i_mobility, fld_norm)
       diffusion_c = LT_get_col(ST_td_tbl, i_diffusion, fld_norm)

       ! The 0.5 is here because of the explicit trapezoidal rule
       dt_vec(ST_ix_cfl) = min(dt_vec(ST_ix_cfl), &
            0.5_dp/sum(abs(fld * mobility) / box%dr))

       ! Dielectric relaxation time
       dt_vec(ST_ix_drt) = min(dt_vec(ST_ix_drt), &
            UC_eps0 / (UC_elem_charge * mobility * &
            max(box%cc(IJK, i_electron), epsilon(1.0_dp))))

       ! Diffusion condition
       dt_vec(ST_ix_diff) = min(dt_vec(ST_ix_diff), &
            0.25_dp * box%dr**2 / max(diffusion_c, epsilon(1.0_dp)))
    end do; CLOSE_DO

  end function get_max_dt

  !> Compute electric field on the tree. First perform multigrid to get electric
  ! potential, then take numerical gradient to geld field.
  subroutine compute_electric_field(tree, n_cycles, have_guess)
    use m_units_constants
    type(a$D_t), intent(inout) :: tree
    integer, intent(in)       :: n_cycles
    logical, intent(in)       :: have_guess
    real(dp), parameter       :: fac = UC_elem_charge / UC_eps0
    integer                   :: lvl, i, id, nc

    nc = tree%n_cell

    ! Set the source term (rhs)
    !$omp parallel private(lvl, i, id)
    do lvl = 1, tree%highest_lvl
       !$omp do
       do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          tree%boxes(id)%cc(DTIMES(:), i_rhs) = fac * (&
               tree%boxes(id)%cc(DTIMES(:), i_electron) - &
               tree%boxes(id)%cc(DTIMES(:), i_pos_ion))
       end do
       !$omp end do nowait
    end do
    !$omp end parallel

    call ST_set_voltage(ST_time)

    ! Perform n_cycles fmg cycles (logicals: store residual, first call)
    do i = 1, n_cycles
       call mg$D_fas_fmg(tree, mg, .false., have_guess .or. i > 1)
    end do

    ! Compute field from potential
    call a$D_loop_box(tree, electric_field_from_potential)

    ! Set the field norm also in ghost cells
    call a$D_gc_tree(tree, i_electric_fld, a$D_gc_interp, a$D_bc_neumann_zero)
  end subroutine compute_electric_field

  !> Compute electric field from electrical potential
  subroutine electric_field_from_potential(box)
    type(box$D_t), intent(inout) :: box
    integer                     :: nc
    real(dp)                    :: inv_dr

    nc     = box%n_cell
    inv_dr = 1 / box%dr

#if $D == 2
    box%fc(1:nc+1, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, i_phi) - box%cc(1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, i_phi) - box%cc(1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 2, electric_fld))**2)
#elif $D == 3
    box%fc(1:nc+1, 1:nc, 1:nc, 1, electric_fld) = inv_dr * &
         (box%cc(0:nc, 1:nc, 1:nc, i_phi) - &
         box%cc(1:nc+1, 1:nc, 1:nc, i_phi))
    box%fc(1:nc, 1:nc+1, 1:nc, 2, electric_fld) = inv_dr * &
         (box%cc(1:nc, 0:nc, 1:nc, i_phi) - &
         box%cc(1:nc, 1:nc+1, 1:nc, i_phi))
    box%fc(1:nc, 1:nc, 1:nc+1, 3, electric_fld) = inv_dr * &
         (box%cc(1:nc, 1:nc, 0:nc, i_phi) - &
         box%cc(1:nc, 1:nc, 1:nc+1, i_phi))

    box%cc(1:nc, 1:nc, 1:nc, i_electric_fld) = 0.5_dp * sqrt(&
         (box%fc(1:nc, 1:nc, 1:nc, 1, electric_fld) + &
         box%fc(2:nc+1, 1:nc, 1:nc, 1, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 2, electric_fld) + &
         box%fc(1:nc, 2:nc+1, 1:nc, 2, electric_fld))**2 + &
         (box%fc(1:nc, 1:nc, 1:nc, 3, electric_fld) + &
         box%fc(1:nc, 1:nc, 2:nc+1, 3, electric_fld))**2)
#endif

  end subroutine electric_field_from_potential

  !> This fills ghost cells near physical boundaries for the potential
  subroutine sides_bc_potential(box, nb, iv, bc_type)
    type(box$D_t), intent(inout) :: box
    integer, intent(in)         :: nb ! Direction for the boundary condition
    integer, intent(in)         :: iv ! Index of variable
    integer, intent(out)        :: bc_type ! Type of boundary condition
    integer                     :: nc

    nc = box%n_cell

    select case (nb)
#if $D == 2
    case (a$D_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, iv) = 0
    case (a$D_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, iv) = 0
    case (a$D_neighb_lowy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc,    0, iv) = 0
    case (a$D_neighb_highy)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, nc+1, iv) = ST_applied_voltage
#elif $D == 3
    case (a3_neighb_lowx)
       bc_type = af_bc_neumann
       box%cc(   0, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_highx)
       bc_type = af_bc_neumann
       box%cc(nc+1, 1:nc, 1:nc, iv) = 0
    case (a3_neighb_lowy)
       bc_type = af_bc_neumann
       box%cc(1:nc,    0, 1:nc, iv) = 0
    case (a3_neighb_highy)
       bc_type = af_bc_neumann
       box%cc(1:nc, nc+1, 1:nc, iv) = 0
    case (a3_neighb_lowz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc,    0, iv) = 0
    case (a3_neighb_highz)
       bc_type = af_bc_dirichlet
       box%cc(1:nc, 1:nc, nc+1, iv) = ST_applied_voltage
#endif
    end select

  end subroutine sides_bc_potential

  !> Compute the electron fluxes due to drift and diffusion
  subroutine fluxes_koren(boxes, id)
    use m_flux_schemes
    type(box$D_t), intent(inout) :: boxes(:)
    integer, intent(in)          :: id
    real(dp)                     :: inv_dr, fld
    real(dp), allocatable        :: v(DTIMES(:), :)
    real(dp), allocatable        :: dc(DTIMES(:), :)
    real(dp), allocatable        :: cc(DTIMES(:))
    integer                      :: nc, n, m
#if $D == 3
    integer                      :: l
#endif

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    allocate(v(DTIMES(1:nc+1), $D))
    allocate(dc(DTIMES(1:nc+1), $D))
    allocate(cc(DTIMES(-1:nc+2)))

    ! Fill ghost cells
    call a$D_gc_box(boxes, id, i_electron, a$D_gc_interp_lim, a$D_bc_neumann_zero)
    call a$D_gc2_box(boxes, id, i_electron, a$D_gc2_prolong_linear, &
         a$D_bc2_neumann_zero, cc, nc)

    ! We use the average field to compute the mobility and diffusion coefficient
    ! at the interface
    do n = 1, nc+1
       do m = 1, nc
#if $D == 2
          fld       = 0.5_dp * (boxes(id)%cc(n-1, m, i_electric_fld) + &
               boxes(id)%cc(n, m, i_electric_fld))
          v(n, m, 1)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
               boxes(id)%fc(n, m, 1, electric_fld)
          dc(n, m, 1) = LT_get_col(ST_td_tbl, i_diffusion, fld)

          fld       = 0.5_dp * (boxes(id)%cc(m, n-1, i_electric_fld) + &
               boxes(id)%cc(m, n, i_electric_fld))
          v(m, n, 2)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
               boxes(id)%fc(m, n, 2, electric_fld)
          dc(m, n, 2) = LT_get_col(ST_td_tbl, i_diffusion, fld)
#elif $D == 3
          do l = 1, nc
             fld = 0.5_dp * (&
                  boxes(id)%cc(n-1, m, l, i_electric_fld) + &
                  boxes(id)%cc(n, m, l, i_electric_fld))
             v(n, m, l, 1)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
                  boxes(id)%fc(n, m, l, 1, electric_fld)
             dc(n, m, l, 1) = LT_get_col(ST_td_tbl, i_diffusion, fld)

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, n-1, l, i_electric_fld) + &
                  boxes(id)%cc(m, n, l, i_electric_fld))
             v(m, n, l, 2)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
                  boxes(id)%fc(m, n, l, 2, electric_fld)
             dc(m, n, l, 2) = LT_get_col(ST_td_tbl, i_diffusion, fld)

             fld = 0.5_dp * (&
                  boxes(id)%cc(m, l, n-1, i_electric_fld) + &
                  boxes(id)%cc(m, l, n, i_electric_fld))
             v(m, l, n, 3)  = -LT_get_col(ST_td_tbl, i_mobility, fld) * &
                  boxes(id)%fc(m, l, n, 3, electric_fld)
             dc(m, l, n, 3) = LT_get_col(ST_td_tbl, i_diffusion, fld)
          end do
#endif
       end do
    end do

    call flux_koren_$Dd(cc, v, nc, 2)
    call flux_diff_$Dd(cc, dc, inv_dr, nc, 2)

    boxes(id)%fc(DTIMES(:), :, flux_elec) = v + dc
  end subroutine fluxes_koren

  !> Take average of new and old electron/ion density for explicit trapezoidal rule
  subroutine average_density(box)
    type(box$D_t), intent(inout) :: box
    box%cc(DTIMES(:), i_electron) = 0.5_dp *  &
         (box%cc(DTIMES(:), i_electron) + box%cc(DTIMES(:), i_electron_old))
    box%cc(DTIMES(:), i_pos_ion) = 0.5_dp * &
         (box%cc(DTIMES(:), i_pos_ion)  + box%cc(DTIMES(:), i_pos_ion_old))
  end subroutine average_density

  !> Advance solution over dt based on the fluxes / source term, using forward Euler
  subroutine update_solution(box, dt)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)         :: dt(:)
    real(dp)                     :: inv_dr, src, fld
    real(dp)                     :: alpha, eta, sflux, mu
#if $D == 2
    real(dp)                     :: rfac(2)
    integer                      :: ioff
#endif
    integer                      :: IJK, nc
    type(LT_loc_t)               :: loc

    nc     = box%n_cell
    inv_dr = 1/box%dr
#if $D == 2
    ioff   = (box%ix(1)-1) * nc
#endif

    do KJI_DO(1,nc)
          fld   = box%cc(IJK, i_electric_fld)
          loc   = LT_get_loc(ST_td_tbl, fld)
          alpha = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          eta   = LT_get_col_at_loc(ST_td_tbl, i_eta, loc)
          mu    = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

          ! Contribution of flux
#if $D == 2
          if (ST_cylindrical) then
             ! Weighting of flux contribution for cylindrical coordinates
             rfac(:) = [i+ioff-1, i+ioff] / (i+ioff-0.5_dp)
          else
             rfac(:) = 1.0_dp
          end if

          sflux = (box%fc(i, j, 2, flux_elec) - box%fc(i, j+1, 2, flux_elec) + &
               rfac(1) * box%fc(i, j, 1, flux_elec) - &
               rfac(2) * box%fc(i+1, j, 1, flux_elec)) * inv_dr * dt(1)
#elif $D == 3
          sflux = (sum(box%fc(i, j, k, :, flux_elec) - &
               box%fc(i+1, j, k, 1, flux_elec) - &
               box%fc(i, j+1, k, 2, flux_elec) - &
               box%fc(i, j, k+1, 3, flux_elec))) * inv_dr * dt(1)
#endif

          ! Source term
          src = fld * mu * box%cc(IJK, i_electron) * (alpha - eta) * dt(1)

          if (ST_photoi_enabled) src = src + box%cc(IJK, i_photo) * dt(1)

          ! Add flux and source term
          box%cc(IJK, i_electron) = box%cc(IJK, i_electron) + sflux + src

          ! Add source term
          box%cc(IJK, i_pos_ion)  = box%cc(IJK, i_pos_ion) + src

    end do; CLOSE_DO
  end subroutine update_solution

  !> Sets the photoionization
  subroutine set_photoionization(tree, eta, num_photons, dt)
    use m_units_constants

    type(a$D_t), intent(inout) :: tree
    real(dp), intent(in)      :: eta
    real(dp), intent(in), optional :: dt
    integer, intent(in)       :: num_photons
    real(dp), parameter       :: p_quench = 30.0e-3_dp
    real(dp)                  :: quench_fac

    ! Compute quench factor, because some excited species will be quenched by
    ! collisions, preventing the emission of a UV photon
    quench_fac = p_quench / (ST_gas_pressure + p_quench)

    ! Set photon production rate per cell, which is proportional to the
    ! ionization rate.
    call a$D_loop_box_arg(tree, set_photoionization_rate, [eta * quench_fac], .true.)

#if $D == 2
    call photoi_set_src_$Dd(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_photo, i_photo, 0.25e-3_dp, .true., ST_cylindrical, 1e-9_dp, dt)
#elif $D == 3
    call photoi_set_src_$Dd(tree, ST_photoi_tbl, ST_rng, num_photons, &
         i_photo, i_photo, 0.25e-3_dp, .true., 1e-9_dp, dt)
#endif

  end subroutine set_photoionization

  !> Sets the photoionization_rate
  subroutine set_photoionization_rate(box, coeff)
    type(box$D_t), intent(inout) :: box
    real(dp), intent(in)        :: coeff(:)
    integer                     :: IJK, nc
    real(dp)                    :: fld, alpha, mobility, tmp
    type(LT_loc_t)              :: loc

    nc = box%n_cell

    do KJI_DO(1,nc)
          fld      = box%cc(IJK, i_electric_fld)
          loc      = LT_get_loc(ST_td_tbl, fld)
          alpha    = LT_get_col_at_loc(ST_td_tbl, i_alpha, loc)
          mobility = LT_get_col_at_loc(ST_td_tbl, i_mobility, loc)

          tmp = fld * mobility * alpha * box%cc(IJK, i_electron) * coeff(1)
          if (tmp < 0) tmp = 0
          box%cc(IJK, i_photo) = tmp
    end do; CLOSE_DO
  end subroutine set_photoionization_rate

  !> For each box that gets refined, set data on its children using this routine
  subroutine prolong_to_new_boxes(tree, ref_info)
    use m_a$D_prolong
    type(a$D_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       !$omp parallel do private(id, p_id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent
          call a$D_prolong_sparse(tree%boxes(p_id), tree%boxes(id), i_electron)
          call a$D_prolong_sparse(tree%boxes(p_id), tree%boxes(id), i_pos_ion)
          call a$D_prolong_sparse(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do
       !$omp end parallel do

       !$omp parallel do private(id)
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          call a$D_gc_box(tree%boxes, id, i_electron, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero)
          call a$D_gc_box(tree%boxes, id, i_pos_ion, &
               a$D_gc_interp_lim, a$D_bc_neumann_zero)
          call a$D_gc_box(tree%boxes, id, i_phi, &
               mg%sides_rb, mg%sides_bc)
       end do
       !$omp end parallel do
     end do
  end subroutine prolong_to_new_boxes

end program streamer_$Dd