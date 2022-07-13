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
  character(len=100) :: fname
  type(af_t)         :: tree ! This contains the full grid information
  type(mg_t)         :: mg   ! Multigrid option struct
  type(ref_info_t)   :: refine_info

  !Indices of cell centered variables
  integer :: i_phi_1 !Phi functions in the improved Eddington Approximation
  integer :: i_phi_2 !Roughly analagous to the electric potential in Poisson's Equation
  integer :: i_elec !Electron density
  integer :: i_elec_old !Old electron density for time-stepping
  integer :: i_photon !photon desnsity
  integer :: i_ph_src !Photoionization source term

  call af_add_cc_variable(tree, "elec", ix=i_elec, n_copies=2)
  i_elec_old = af_find_cc_variable(tree, "elec_2")
  call af_add_cc_variable(tree, "phi_1", ix=i_phi_1)
  call af_add_cc_variable(tree, "phi_2", ix=i_phi_2)
  call af_add_cc_variable(tree, "ph_src", ix=i_ph_src)
  call af_add_cc_variable(tree, "photon", ix=i_photon)

  ! Simulation parameters
  real(dp), parameter :: end_time      = 10e-9_dp
  real(dp), parameter :: dt_output     = 20e-11_dp
  real(dp), parameter :: dt_max        = 20e-11_dp
  integer, parameter  :: ref_per_steps = 2
  integer, parameter  :: box_size      = 8

  !photo-ionization parameters
  real(dp), parameter :: ph_quenching = 1
  real(dp), parameter :: ph_chi_min = 0.035 !/(torr*cm)
  real(dp), parameter :: ph_chi_max = 2 !/(torr*cm)
  real(dp), parameter :: gamma_1 = 5/7 * (1 - 3 * sqrt(6 / 5))
  real(dp), parameter :: gamma_1 = 5/7 * (1 + 3 * sqrt(6 / 5))
  !Fit parameters
  !Reported in A Bourdon et al 2007 Plasma Sources Sci. Technol 16 656
  real(dp), parameter :: A_j = (/ 0.0067_dp,0.0346_dp,0.3059_dp /)
  real(dp), parameter :: lambda_j = (/ 0.0447_dp,0.1121_dp,0.5994_dp /)

  ! Computational domain
  real(dp), parameter :: domain_length = 2e-3_dp
  real(dp), parameter :: refine_max_dx = 1e-3_dp
  real(dp), parameter :: refine_min_dx = 1e-9_dp

  !Set multigrid options


contains
  subroutine compute_ph_src

  end subroutine compute_ph_src

  !Compute the photon desnsity from the functions phi_1 and phi_2
  !See pp. 662 of Bourdon et al 2007
  subroutine photon_from_phi(box)
    type(box), intent(inout) :: box 
    integer                  :: nc

    nc = box%n_cell

    box%cc(1:nc, 1:nc, i_photon) = (gamma_2 * box%cc(1:nc, 1:nc, i_phi_1) - &
      gamma_1 * box%cc(1:nc, 1:nc, i_phi__2)) / (gamma_2 - gamma_1)
  end subroutine ph_src_from_phi



end program photoionization
