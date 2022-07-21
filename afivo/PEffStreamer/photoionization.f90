#include "../src/cpp_macros.h"

program photoionization
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
  type(af_t)         :: tree ! This contains the full grid information
  type(mg_t), dimension(3)         :: mg   ! Multigrid option struct for each Helmholtz group
  type(ref_info_t)   :: refine_info

  !Indices of cell centered variables
  integer :: i_ph_src !Photoionization source term
  integer :: i_ph_dist !Distribution of photoionizing photons

  integer :: i_psi_1 !Isotropic component of photon distribution function
  integer :: i_psi_2 !For the jth wavelength group
  integer :: i_psi_3 !j = 1, 2, 3

  integer :: i_rhs_1 !Right hand variable for Eddington Approximation
  integer :: i_rhs_2
  integer :: i_rhs_3

  integer :: i_tmp !Temporary variable

  integer :: psi_indices(3) !arrays to hold indices of psi for easy access in loops
  integer :: rhs_indices(3)



  ! Simulation parameters
  real(dp), parameter :: end_time      = 10e-9_dp
  real(dp), parameter :: dt_output     = 20e-11_dp
  real(dp), parameter :: dt_max        = 20e-11_dp
  integer, parameter  :: ref_per_steps = 2
  integer, parameter  :: box_size      = 8
  real(dp), parameter :: pi            = 3.14159265359



  !photo-ionization parameters
  real(dp), parameter :: ph_quenching = 1
  real(dp), parameter :: ph_chi_min = 0.035 !/(torr*cm)
  real(dp), parameter :: ph_chi_max = 2 !/(torr*cm)
  real(dp), parameter :: I_0 = 3.5e22_dp !/(cm**3*s)
  real(dp), parameter :: p_o2 = 150 !Partial pressure of oxygen; 150 at ATP
  real(dp), parameter :: xi = 0.1 !photoionization efficiency
  !Fit parameters
  !Reported in A Bourdon et al 2007 Plasma Sources Sci. Technol 16 656
  real(dp), parameter :: A(3) = (/ 0.0067_dp,0.0346_dp,0.3059_dp /)
  real(dp), parameter :: lambda(3) = (/ 0.0447_dp,0.1121_dp,0.5994_dp /)

  ! Computational domain
  real(dp), parameter :: domain_length = 2e-3_dp
  real(dp), parameter :: refine_max_dx = 1e-3_dp
  real(dp), parameter :: refine_min_dx = 1e-9_dp
  real(dp), parameter :: sigma         =  0.01 !cm
  real(dp), parameter :: x_0           = 0.5 * domain_length
  real(dp), parameter :: y_0           = 0.5 * domain_length

  call af_add_cc_variable(tree, "ph_src", ix=i_ph_src)
  call af_add_cc_variable(tree, "ph_dist", ix=i_ph_dist)
  call af_add_cc_variable(tree, "psi_1", ix=i_psi_1)
  call af_add_cc_variable(tree, "psi_2", ix=i_psi_2)
  call af_add_cc_variable(tree, "psi_3", ix=i_psi_3)
  call af_add_cc_variable(tree, "rhs_1", ix=i_rhs_1)
  call af_add_cc_variable(tree, "rhs_2", ix=i_rhs_2)
  call af_add_cc_variable(tree, "rhs_3", ix=i_rhs_3)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)

  psi_indices(1) = i_psi_1
  psi_indices(2) = i_psi_2
  psi_indices(3) = i_psi_3

  rhs_indices(1) = i_rhs_1
  rhs_indices(2) = i_rhs_2
  rhs_indices(3) = i_rhs_3

  call af_init(tree, &
               box_size, &
               [domain_length, domain_length], &
               [box_size, box_size], &
               coord=af_cyl)

  !Set multigrid options
  mg(1)%i_phi     = i_psi_1
  mg(1)%i_rhs     = i_rhs_1
  mg(1)%i_tmp     = i_tmp
  mg(1)%helmholtz_lambda = (3 * (lambda(1) * p_o2) ** 2)

  mg(2)%i_phi     = i_psi_2
  mg(2)%i_rhs     = i_rhs_2
  mg(2)%i_tmp     = i_tmp
  mg(2)%helmholtz_lambda = (3 * (lambda(2) * p_o2) ** 2)

  mg(3)%i_phi     = i_psi_3
  mg(3)%i_rhs     = i_rhs_3
  mg(3)%i_tmp     = i_tmp
  mg(3)%helmholtz_lambda = (3 * (lambda(3) * p_o2) ** 2)


  do i=1, 3
    if (i == 1) then
      mg(i)%sides_bc => af_bc_dirichlet_zero
    else
      mg(i)%sides_bc => af_bc_dirichlet_zero
    end if

    call mg_init(tree, mg(i))
  end do

  !Set initial conditions
  do i=1, 5
    print *, "Initializing", i
    call af_loop_box(tree, set_ph_dist)
    call compute_ph_src(tree, .false.)
    call af_adjust_refinement(tree, refinement_routine, refine_info)


    if (refine_info%n_add == 0 .and. refine_info%n_rm == 0) then
      print *, "Nothing to add: exiting loop prematuraly"
      exit
    end if
  end do

  call af_loop_box(tree, set_ph_dist)
  call compute_ph_src(tree, .true.)

  !write to silo
  call af_write_silo(tree, "photoionization_init")

contains
  !Refinement routine
  !Just going to use the one from the Poisson_Helmholtz example and see how it works
  !If its shit, I will try to think something better up
  subroutine refinement_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(box%n_cell, box%n_cell)
    integer                 :: i, j, nc
    real(dp)                :: xy(2), dr2, drhs

    nc = box%n_cell
    dr2 = maxval(box%dr)**2
    do i=1, nc
      do j=1, nc
        xy = af_r_cc(box, [i,j])
        drhs = dr2 * max(box%cc(i, j, i_rhs_1), box%cc(i, j, i_rhs_2), box%cc(i, j, i_rhs_3))
        if (abs(drhs) > 1e-4 .and. box%lvl<5) then
          cell_flags(i,j) = af_do_ref
        else
          cell_flags(i,j) = af_keep_ref
        end if
      end do
    end do
  end subroutine refinement_routine

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
      call mg_fas_fmg(tree, mg(j), .false., have_guess)
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

  subroutine set_ph_dist(box)
    type(box_t), intent(inout) :: box
    real(dp)                   :: xy(2)
    integer                    :: i, j, nc

    !print *, "set photon distribution"
    nc = box%n_cell

    do j=1, nc
      do i=1, nc
        xy = af_r_cc(box, [i,j])
        !print *,  I_0 * exp(-((xy(1) - x_0)**2 + (xy(2) - y_0)**2) / sigma**2)
        box%cc(i, j, i_ph_dist) = I_0 * exp(-((xy(1))**2 + (xy(2) - y_0)**2) / sigma**2)
      end do
    end do
  end subroutine set_ph_dist

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
    real(dp)                :: s_ph, xy_1(2), xy_2(2), R, vol, dr(2)

    !Dirichlet boundary conditions
    bc_type = af_bc_dirichlet
    !print *, "Setting boundary conditions"


    do j = 1, box%n_cell
      xy_1(1) = coords(1, j)
      xy_1(2) = coords(2, j)
      !print *, xy_1
      s_ph = 0

      do lvl = 1, tree%highest_lvl
        do i = 1, size(tree%lvls(lvl)%leaves)
          id = tree%lvls(lvl)%leaves(i)
          nc = tree%boxes(id)%n_cell
          !print *, "Level, i: ", lvl, i

          !Solve for boundary conditions using Zheleznyak integral
          do l = 1, nc
            do m = 1, nc

              xy_2 = af_r_cc(tree%boxes(id), [l, m])
              R = sqrt((xy_2(1) - xy_1(1))**2 + (xy_2(2) - xy_1(2))**2)
              !print *, R, xy_1
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

  real(dp) function f(R)
    real(dp), intent(in) :: R

    f = (exp(-ph_chi_min * p_o2 * R) - exp(-ph_chi_max * p_o2 * R)) / (R * log(ph_chi_max / ph_chi_min))
  end function f
end program photoionization
