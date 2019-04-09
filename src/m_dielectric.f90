#include "../afivo/src/cpp_macros.h"
!> Module with settings and routines to handle dielectrics
!>
!> @todo Make sure species densities are initially zero inside the dielectric
!> @todo Use special prolongation for multigrid when there are surface charges
module m_dielectric
  use m_af_all

  implicit none
  private

  type surface_data_t
     logical :: in_use = .false.
     integer :: id_gas
     integer :: id_diel
     integer :: direction
     real(dp) :: eps
#if NDIM == 2
     !> Charge density on the dielectric
     real(dp), allocatable :: charge(:, :)
     !> Photon flux on the dielectric
     real(dp), allocatable :: photons(:)
#elif NDIM == 3
     !> Charge density on the dielectric
     real(dp), allocatable :: charge(:, :, :)
     !> Photon flux on the dielectric
     real(dp), allocatable :: photons(:, :)
#endif
  end type surface_data_t

  integer                           :: num_surfaces = 0
  type(surface_data_t), allocatable, public :: surface_list(:)
  integer                           :: n_removed_surfaces = 0
  integer, allocatable              :: removed_surfaces(:)
  integer, allocatable, public, protected :: box_id_to_surface_id(:)

  ! Maximum travel distance for testing boundary intersection
  real(dp), protected :: photon_step_length = 1.0e-3_dp

  public :: dielectric_initialize
  public :: dielectric_allocate
  public :: dielectric_get_surfaces
  public :: dielectric_update_after_refinement
  public :: dielectric_fix_refine
  public :: dielectric_rearrange_charge
  public :: dielectric_adjust_field
  public :: dielectric_update_surface_charge
  public :: dielectric_photon_emission
  public :: dielectric_combine_substeps
  public :: dielectric_photon_absorption

contains

  subroutine dielectric_initialize(tree, cfg)
    use m_config
    type(af_t), intent(in) :: tree
    type(CFG_t), intent(inout) :: cfg

    call CFG_add_get(cfg, "dielectric%photon_step_length", &
         photon_step_length, &
         "Maximum travel distance for testing boundary intersection")
  end subroutine dielectric_initialize

  subroutine dielectric_allocate(tree)
    use m_config
    type(af_t), intent(in) :: tree

    allocate(box_id_to_surface_id(tree%box_limit))

    ! Maximum number of surfaces
    allocate(surface_list(tree%box_limit/2))
    allocate(removed_surfaces(tree%box_limit/2))

    box_id_to_surface_id(:) = -1
  end subroutine dielectric_allocate

  integer function get_new_surface()
    if (n_removed_surfaces > 0) then
       get_new_surface = removed_surfaces(n_removed_surfaces)
       n_removed_surfaces = n_removed_surfaces - 1
    else
       num_surfaces = num_surfaces + 1
       get_new_surface = num_surfaces
    end if
  end function get_new_surface

  subroutine dielectric_get_surfaces(tree)
    use m_advance_base
    type(af_t), intent(in) :: tree
    integer                :: lvl, i, id, nb, nc, ix
    integer                :: id_diel
    real(dp)               :: eps

    nc = tree%n_cell

    ! Locate all boxes at the boundary of the dielectric
    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb, eps)

          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%in_use = .true.
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb
             surface_list(ix)%eps = eps

             if (.not. allocated(surface_list(ix)%charge)) then
#if NDIM == 2
                allocate(surface_list(ix)%charge(nc, advance_num_states))
                allocate(surface_list(ix)%photons(nc))
#elif NDIM == 3
                allocate(surface_list(ix)%charge(nc, nc, advance_num_states))
                allocate(surface_list(ix)%photons(nc, nc))
#endif
             end if
             surface_list(ix)%charge = 0.0_dp
             surface_list(ix)%photons = 0.0_dp

          end if
       end do
    end do

  end subroutine dielectric_get_surfaces

  subroutine find_dielectric_neighbor(tree, id, id_diel, nb, nb_eps)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: id
    integer, intent(out)   :: id_diel
    integer, intent(out)   :: nb
    real(dp), intent(out)  :: nb_eps
    integer                :: nb_id
    real(dp), parameter    :: eps_threshold = 1e-8_dp

    id_diel = -1
    nb      = -1
    nb_eps  = -1.0_dp

    if (eps_changes_in_box(tree%boxes(id)) .and. &
         tree%boxes(id)%cc(DTIMES(1), i_eps) < 1 + eps_threshold) then

       do nb = 1, af_num_neighbors
          nb_id = tree%boxes(id)%neighbors(nb)
          if (nb_id <= af_no_box) cycle

          nb_eps = tree%boxes(nb_id)%cc(DTIMES(1), i_eps)
          if (nb_eps > 1 + eps_threshold) then
             id_diel = nb_id
             exit
          end if
       end do
    end if

  end subroutine find_dielectric_neighbor

  logical function eps_changes_in_box(box)
    use m_streamer
    type(box_t), intent(in) :: box
    real(dp)                :: a, b
#if NDIM == 2
    a = minval(box%cc(:, :, i_eps))
    b = maxval(box%cc(:, :, i_eps))
#elif NDIM == 3
    a = minval(box%cc(:, :, :, i_eps))
    b = maxval(box%cc(:, :, :, i_eps))
#endif

    eps_changes_in_box = (b > a)
  end function eps_changes_in_box

  subroutine dielectric_update_after_refinement(tree, ref_info)
    use m_advance_base
    type(af_t), intent(in)       :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, id_diel, nb, ix
    integer                      :: nc
    real(dp)                     :: eps

    nc = tree%n_cell

    ! Handle removed surfaces
    do i = 1, size(ref_info%rm)
       id = ref_info%rm(i)
       ix = box_id_to_surface_id(id)
       if (ix > 0) then
          call restrict_surface_to_parent(tree, surface_list(ix))

          n_removed_surfaces = n_removed_surfaces + 1
          removed_surfaces(n_removed_surfaces) = ix
          box_id_to_surface_id(surface_list(ix)%id_gas) = -1
          box_id_to_surface_id(surface_list(ix)%id_diel) = -1
          surface_list(ix)%in_use = .false.
       end if
    end do

    ! Add new surfaces
    do lvl = 1, tree%highest_lvl
       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)

          call find_dielectric_neighbor(tree, id, id_diel, nb, eps)

          if (id_diel > 0) then
             ix = get_new_surface()
             surface_list(ix)%in_use = .true.
             surface_list(ix)%id_gas = id
             surface_list(ix)%id_diel = id_diel
             surface_list(ix)%direction = nb
             surface_list(ix)%eps = eps

             box_id_to_surface_id(id) = ix
             box_id_to_surface_id(id_diel) = ix

             if (.not. allocated(surface_list(ix)%charge)) then
#if NDIM == 2
                allocate(surface_list(ix)%charge(nc, advance_num_states))
                allocate(surface_list(ix)%photons(nc))
#elif NDIM == 3
                allocate(surface_list(ix)%charge(nc, nc, advance_num_states))
                allocate(surface_list(ix)%photons(nc, nc))
#endif
             end if
             surface_list(ix)%charge = 0.0_dp
             surface_list(ix)%photons = 0.0_dp

             call prolong_surface_from_parent(tree, surface_list(ix))
          end if
       end do
    end do

  end subroutine dielectric_update_after_refinement

  subroutine prolong_surface_from_parent(tree, surface)
    type(af_t), intent(in)              :: tree
    type(surface_data_t), intent(inout) :: surface
    integer                             :: id, p_id, ix_p, ix_offset(NDIM)
    integer                             :: nc, dix

    ! Find parent surface
    id = surface%id_gas
    p_id = tree%boxes(id)%parent
    ix_p = box_id_to_surface_id(p_id)

    if (ix_p > 0) then
       nc        = tree%n_cell
       ix_offset = af_get_child_offset(tree%boxes(id))

       ! Determine dimension along surface
#if NDIM == 2
       if (af_neighb_dim(surface%direction) == 1) then
          dix = ix_offset(2)
       else
          dix = ix_offset(1)
       end if

       ! Simply copy the values from the parent
       surface%charge(1:nc:2, 1) = &
            surface_list(ix_p)%charge(dix+1:dix+nc/2, 1)
       surface%charge(2:nc:2, 1) = surface%charge(1:nc:2, 1)

       surface%photons(1:nc:2) = &
            surface_list(ix_p)%photons(dix+1:dix+nc/2)
       surface%photons(2:nc:2) = surface%photons(1:nc:2)
#elif NDIM == 3
       error stop
#endif
    end if

  end subroutine prolong_surface_from_parent

  subroutine restrict_surface_to_parent(tree, surface)
    type(af_t), intent(in)           :: tree
    type(surface_data_t), intent(in) :: surface
    integer                          :: id, p_id, ix_p, ix_offset(NDIM)
    integer                          :: nc, dix

    ! Find parent surface
    id   = surface%id_gas
    p_id = tree%boxes(id)%parent
    ix_p = box_id_to_surface_id(p_id)

    if (ix_p > 0) then
       nc        = tree%n_cell
       ix_offset = af_get_child_offset(tree%boxes(id))

       ! Determine dimension along surface
#if NDIM == 2
       if (af_neighb_dim(surface%direction) == 1) then
          dix = ix_offset(2)
       else
          dix = ix_offset(1)
       end if

       ! Average the value on the children
       surface_list(ix_p)%charge(dix+1:dix+nc/2, 1) = 0.5_dp * ( &
            surface%charge(1:nc:2, 1) + surface%charge(2:nc:2, 1))
       surface_list(ix_p)%photons(dix+1:dix+nc/2) = 0.5_dp * ( &
            surface%photons(1:nc:2) + surface%photons(2:nc:2))
#elif NDIM == 3
       error stop
#endif
    end if

  end subroutine restrict_surface_to_parent

  subroutine dielectric_fix_refine(tree, ref_flags)
    use m_streamer
    type(af_t), intent(in) :: tree
    integer, intent(inout) :: ref_flags(:)
    integer                :: lvl, i, id, nb_id, ix

    do lvl = 1, tree%highest_lvl
       do i = 1, size(tree%lvls(lvl)%ids)
          id = tree%lvls(lvl)%ids(i)

          ix = box_id_to_surface_id(id)
          if (ix > 0) then
             if (id == surface_list(ix)%id_gas) then
                nb_id = surface_list(ix)%id_diel
                ref_flags([id, nb_id]) = maxval(ref_flags([id, nb_id]))
             end if
          end if
       end do
    end do
  end subroutine dielectric_fix_refine

  subroutine dielectric_rearrange_charge(tree, s_in)
    type(af_t), intent(inout) :: tree
    integer, intent(in)       :: s_in !< Time state to use
    integer                   :: n, id_gas

    do n = 1, num_surfaces
       if (.not. surface_list(n)%in_use) cycle
       id_gas = surface_list(n)%id_gas

       if (af_has_children(tree%boxes(id_gas))) cycle

       call surface_charge_to_rhs(tree%boxes, &
            surface_list(n), tree%n_cell, s_in)
    end do
  end subroutine dielectric_rearrange_charge

  subroutine surface_charge_to_rhs(boxes, surface, nc, s_in)
    use m_streamer
    use m_units_constants
    type(box_t), intent(inout)          :: boxes(:)
    type(surface_data_t), intent(inout) :: surface
    integer, intent(in)                 :: nc
    integer, intent(in)                 :: s_in !< Time state to use
    integer                             :: nb, id_gas, id_diel
    integer                             :: glo(NDIM), ghi(NDIM)
    integer                             :: dlo(NDIM), dhi(NDIM)
    real(dp)                            :: frac_gas, dr, fac

    nb      = surface%direction
    id_gas  = surface%id_gas
    id_diel = surface%id_diel
    dr      = boxes(id_gas)%dr(af_neighb_dim(nb))

    ! Factor to convert surface charge density to rhs, which is defined as
    ! -charge_density / eps0
    fac     = -1 / (dr * UC_eps0)

    ! Get index range for dielectric box
    call af_get_index_bc_inside(af_neighb_rev(nb), nc, dlo, dhi)

    ! Get index range for gas-phase box
    call af_get_index_bc_inside(nb, nc, glo, ghi)

    ! How much of the rhs to put on the gas and dielectric side
    frac_gas  = 1.0_dp / (1.0_dp + surface%eps)

#if NDIM == 2
    boxes(id_gas)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) = &
         boxes(id_gas)%cc(glo(1):ghi(1), glo(2):ghi(2), i_rhs) &
         + frac_gas * fac * &
         reshape(surface%charge(:, 1+s_in), [ghi(1)-glo(1)+1, ghi(2)-glo(2)+1])

    boxes(id_diel)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_rhs) = &
         boxes(id_diel)%cc(dlo(1):dhi(1), dlo(2):dhi(2), i_rhs) &
         + (1-frac_gas) * fac * &
         reshape(surface%charge(:, 1+s_in), [dhi(1)-dlo(1)+1, dhi(2)-dlo(2)+1])
#elif NDIM == 3
    error stop
#endif

  end subroutine surface_charge_to_rhs

  subroutine dielectric_adjust_field(boxes, id, s_in)
    use m_units_constants
    use m_streamer
    type(box_t), intent(inout) :: boxes(:)
    integer, intent(in)        :: id
    integer, intent(in)        :: s_in
    integer                    :: nc, nb, ix
    real(dp)                   :: eps, fac_fld, fac_charge

    nc = boxes(id)%n_cell

    ix = box_id_to_surface_id(id)

    if (ix <= 0) return

    eps = surface_list(ix)%eps

    if (id == surface_list(ix)%id_gas) then
       nb = surface_list(ix)%direction
       fac_fld = 2 * eps / (1 + eps)
       ! This factor is 1/eps0 * 1/(eps_gas + eps_diel)
       fac_charge = 1 / (UC_eps0 * (1 + eps))
    else
       nb = af_neighb_rev(surface_list(ix)%direction)
       fac_fld = 2 / (1 + eps)
       fac_charge = 1 / (UC_eps0 * (1 + eps))
    end if

#if NDIM == 2
    ! Compute fields at the boundaries of the box, where eps can change
    select case (nb)
    case (af_neighb_lowx)
       boxes(id)%fc(1, 1:nc, 1, electric_fld) = fac_fld * &
            boxes(id)%fc(1, 1:nc, 1, electric_fld) + &
            fac_charge * surface_list(ix)%charge(:, 1+s_in)
    case (af_neighb_highx)
       boxes(id)%fc(nc+1, 1:nc, 1, electric_fld) = fac_fld * &
            boxes(id)%fc(nc+1, 1:nc, 1, electric_fld) - &
            fac_charge * surface_list(ix)%charge(:, 1+s_in)
    case (af_neighb_lowy)
       boxes(id)%fc(1:nc, 1, 2, electric_fld) = fac_fld * &
            boxes(id)%fc(1:nc, 1, 2, electric_fld) + &
            fac_charge * surface_list(ix)%charge(:, 1+s_in)
    case (af_neighb_highy)
       boxes(id)%fc(1:nc, nc+1, 2, electric_fld) = fac_fld * &
            boxes(id)%fc(1:nc, nc+1, 2, electric_fld) - &
            fac_charge * surface_list(ix)%charge(:, 1+s_in)
    end select
#elif NDIM == 3
    error stop
    ! box%fc(1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
    !      (box%cc(0, 1:nc, 1:nc, i_phi) - box%cc(1, 1:nc, 1:nc, i_phi)) * &
    !      box%cc(0, 1:nc, 1:nc, i_eps) / &
    !      (box%cc(1, 1:nc, 1:nc, i_eps) + box%cc(0, 1:nc, 1:nc, i_eps))
    ! box%fc(nc+1, 1:nc, 1:nc, 1, electric_fld) = 2 * inv_dr(1) * &
    !      (box%cc(nc, 1:nc, 1:nc, i_phi) - box%cc(nc+1, 1:nc, 1:nc, i_phi)) * &
    !      box%cc(nc+1, 1:nc, 1:nc, i_eps) / &
    !      (box%cc(nc+1, 1:nc, 1:nc, i_eps) + box%cc(nc, 1:nc, 1:nc, i_eps))
    ! box%fc(1:nc, 1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
    !      (box%cc(1:nc, 0, 1:nc, i_phi) - box%cc(1:nc, 1, 1:nc, i_phi)) * &
    !      box%cc(1:nc, 0, 1:nc, i_eps) / &
    !      (box%cc(1:nc, 1, 1:nc, i_eps) + box%cc(1:nc, 0, 1:nc, i_eps))
    ! box%fc(1:nc, nc+1, 1:nc, 2, electric_fld) = 2 * inv_dr(2) * &
    !      (box%cc(1:nc, nc, 1:nc, i_phi) - box%cc(1:nc, nc+1, 1:nc, i_phi)) * &
    !      box%cc(1:nc, nc+1, 1:nc, i_eps) / &
    !      (box%cc(1:nc, nc+1, 1:nc, i_eps) + box%cc(1:nc, nc, 1:nc, i_eps))
    ! box%fc(1:nc, 1:nc, 1, 3, electric_fld) = 2 * inv_dr(3) * &
    !      (box%cc(1:nc, 1:nc, 0, i_phi) - box%cc(1:nc, 1:nc, 1, i_phi)) * &
    !      box%cc(1:nc, 1:nc, 0, i_eps) / &
    !      (box%cc(1:nc, 1:nc, 1, i_eps) + box%cc(1:nc, 1:nc, 0, i_eps))
    ! box%fc(1:nc, 1:nc, nc+1, 3, electric_fld) = 2 * inv_dr(3) * &
    !      (box%cc(1:nc, 1:nc, nc, i_phi) - box%cc(1:nc, 1:nc, nc+1, i_phi)) * &
    !      box%cc(1:nc, 1:nc, nc+1, i_eps) / &
    !      (box%cc(1:nc, 1:nc, nc+1, i_eps) + box%cc(1:nc, 1:nc, nc, i_eps))
#endif

  end subroutine dielectric_adjust_field

  subroutine dielectric_update_surface_charge(box, surface, dt, s_in, s_out)
    use m_units_constants
    use m_streamer
    type(box_t), intent(inout)          :: box
    type(surface_data_t), intent(inout) :: surface
    real(dp), intent(in)                :: dt    !< Time step
    integer, intent(in)                 :: s_in  !< Input state
    integer, intent(in)                 :: s_out !< Output state
    integer                             :: nc
    real(dp)                            :: fac

    nc  = box%n_cell
    fac = dt * UC_elec_charge

    select case (surface%direction)
#if NDIM == 2
    case (af_neighb_lowx)
       surface%charge(:, 1+s_out) = surface%charge(:, 1+s_in) - &
            fac * box%fc(1, 1:nc, 1, flux_elec)
    case (af_neighb_highx)
       surface%charge(:, 1+s_out) = surface%charge(:, 1+s_in) + &
            fac * box%fc(nc+1, 1:nc, 1, flux_elec)
    case (af_neighb_lowy)
       surface%charge(:, 1+s_out) = surface%charge(:, 1+s_in) - &
            fac * box%fc(1:nc, 1, 2, flux_elec)
    case (af_neighb_highy)
       surface%charge(:, 1+s_out) = surface%charge(:, 1+s_in) + &
            fac * box%fc(1:nc, nc+1, 2, flux_elec)
#elif NDIM == 3
    case default
       error stop
#endif
    end select
  end subroutine dielectric_update_surface_charge

  subroutine dielectric_photon_emission(box, surface, dt, s_in, s_out)
    use m_units_constants
    use m_streamer
    type(box_t), intent(inout)          :: box
    type(surface_data_t), intent(inout) :: surface
    real(dp), intent(in)                :: dt    !< Time step
    integer, intent(in)                 :: s_in  !< Input state
    integer, intent(in)                 :: s_out !< Output state
    integer                             :: nc
    real(dp)                            :: fac, dr
    real(dp), parameter                 :: gamma = 0.1_dp

    nc  = box%n_cell
    dr  = box%dr(af_neighb_dim(surface%direction))
    fac = gamma * dt

    select case (surface%direction)
#if NDIM == 2
    case (af_neighb_lowx)
       where (box%fc(1, 1:nc, 1, electric_fld) < 0.0_dp)
          box%cc(1, 1:nc, i_electron+s_out) = &
               box%cc(1, 1:nc, i_electron+s_out) + &
               surface%photons * fac / dr
          surface%charge(:, 1+s_out) = surface%charge(:, 1+s_out) + &
               surface%photons * fac * UC_elem_charge
       end where
    case (af_neighb_highx)
       where (box%fc(nc, 1:nc, 1, electric_fld) > 0.0_dp)
          box%cc(nc, 1:nc, i_electron+s_out) = &
               box%cc(nc, 1:nc, i_electron+s_out) + &
               surface%photons * fac / dr
          surface%charge(:, 1+s_out) = surface%charge(:, 1+s_out) + &
               surface%photons * fac * UC_elem_charge
       end where
    case (af_neighb_lowy)
       where (box%fc(1:nc, 1, 2, electric_fld) < 0.0_dp)
          box%cc(1:nc, 1, i_electron+s_out) = &
               box%cc(1:nc, 1, i_electron+s_out) + &
               surface%photons * fac / dr
          surface%charge(:, 1+s_out) = surface%charge(:, 1+s_out) + &
               surface%photons * fac * UC_elem_charge
       end where
    case (af_neighb_highy)
       where (box%fc(1:nc, nc, 2, electric_fld) > 0.0_dp)
          box%cc(1:nc, nc, i_electron+s_out) = &
               box%cc(1:nc, nc, i_electron+s_out) + &
               surface%photons * fac / dr
          surface%charge(:, 1+s_out) = surface%charge(:, 1+s_out) + &
               surface%photons * fac * UC_elem_charge
       end where
#elif NDIM == 3
    case default
       error stop
#endif
    end select

  end subroutine dielectric_photon_emission

  subroutine dielectric_combine_substeps(tree, in_steps, coeffs, out_step)
    type(af_t), intent(in) :: tree
    integer, intent(in)    :: in_steps(:)
    real(dp), intent(in)   :: coeffs(:)
    integer, intent(in)    :: out_step
    integer                :: n, k
#if NDIM == 2
    real(dp)               :: tmp(tree%n_cell)
#elif NDIM == 3
    real(dp)               :: tmp(tree%n_cell, tree%n_cell)
#endif

    do n = 1, num_surfaces
       if (.not. surface_list(n)%in_use) cycle

       tmp = 0.0_dp
       do k = 1, size(in_steps)
#if NDIM == 2
          tmp = tmp + coeffs(k) * surface_list(n)%charge(:, 1+in_steps(k))
#elif NDIM == 3
          error stop
#endif
       end do
       surface_list(n)%charge(:, 1+out_step) = tmp
    end do
  end subroutine dielectric_combine_substeps

  !> Determine whether and where photons crossed a dielectric, and change their
  !> absorption location to the first cell inside the surface
  subroutine dielectric_photon_absorption(tree, xyz_start, xyz_end, n_dim, &
       n_photons, photon_weight)
    use m_af_types
    use m_af_interp
    use m_streamer
    type(af_t), intent(in)  :: tree
    integer, intent(in)     :: n_dim
    integer, intent(in)     :: n_photons
    real(dp), intent(in)    :: xyz_start(3, n_photons)
    real(dp), intent(inout) :: xyz_end(3, n_photons)
    real(dp), intent(in)    :: photon_weight
    real(dp)                :: xyz(n_dim), dvec(n_dim)
    real(dp)                :: xyz_gas(n_dim), xyz_nogas(n_dim)
    real(dp)                :: xyz_middle(n_dim), eps(1)
    real(dp)                :: travel_distance
    integer                 :: n, n_steps, i, k, n_bisect
    logical                 :: success

    ! Determine the number of bisection steps to find the first cell inside the
    ! dielectric, given a photon step length <= photon_step_length
    n_bisect = -ceiling(&
         log(af_min_dr(tree)/photon_step_length) / log(2.0_dp))

    do n = 1, n_photons
       xyz             = xyz_start(1:n_dim, n)
       dvec            = xyz_end(1:n_dim, n) - xyz_start(1:n_dim, n)
       travel_distance = norm2(dvec)
       n_steps         = ceiling(travel_distance/photon_step_length)

       ! Normalize direction vector to right length. Possible TODO: near the
       ! boundary of the domain, a photon can fly out crossing on a small
       ! piece of a dielectric.
       dvec = dvec / n_steps

       do i = 1, n_steps
          xyz = xyz + dvec

          ! If outside, stop
          if (any(xyz < ST_domain_origin .or. &
               xyz > ST_domain_origin + ST_domain_len)) then
             exit
          end if

          ! Get dielectric permittivity
          eps = af_interp0(tree, xyz, [i_eps], success)
          if (.not. success) error stop "photon unexpectedly outside domain"

          ! If epsilon changes, start doing a search for the intersection point
          if (eps(1) > 1.0_dp) then
             ! Perform bisection to locate first cell inside dielectric
             xyz_gas = xyz - dvec
             xyz_nogas = xyz

             do k = 1, n_bisect
                xyz_middle = 0.5_dp * (xyz_gas + xyz_nogas)
                eps = af_interp0(tree, xyz_middle, [i_eps], success)
                if (.not. success) error stop "photon unexpectedly outside domain"

                if (eps(1) > 1.0_dp) then
                   xyz_nogas = xyz_middle
                else
                   xyz_gas = xyz_middle
                end if
             end do

             call add_to_surface_photons(tree, xyz_nogas, photon_weight)
             xyz_end(1:n_dim, n) = -1e50_dp
             exit
          end if
       end do
    end do
  end subroutine dielectric_photon_absorption

  subroutine add_to_surface_photons(tree, xyz, w)
    type(af_t), intent(in) :: tree
    real(dp), intent(in)   :: xyz(NDIM)
    real(dp), intent(in)   :: w
    type(af_loc_t)         :: loc
    integer                :: ix, id_gas, normal_dim, i
    real(dp)               :: area

    loc = af_get_loc(tree, xyz)

    if (loc%id == -1) error stop "wrong photon location"

    ix = box_id_to_surface_id(loc%id)
    if (ix == -1) error stop "no surface here"

    id_gas = surface_list(ix)%id_gas
    normal_dim = af_neighb_dim(surface_list(ix)%direction)

#if NDIM == 2
    if (normal_dim == 1) then
       i = loc%ix(2)
       area = tree%boxes(id_gas)%dr(2)
    else
       i = loc%ix(1)
       area = tree%boxes(id_gas)%dr(1)
    end if

    surface_list(ix)%photons(i) = surface_list(ix)%photons(i) + w/area
#elif NDIM == 3
    error stop
#endif
  end subroutine add_to_surface_photons

end module m_dielectric
