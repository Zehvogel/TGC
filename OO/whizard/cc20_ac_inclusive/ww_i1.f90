! File generated automatically by O'Mega 3.1.4 release Nov 08 2023
!
!   /cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/whizard/3.1.4-zjmc7r/bin/omega_SM_ac.opt -o ww_i1.f90 -target:whizard -target:parameter_module parameters_SM_ac -target:module opr_ww_i1 -target:md5sum 84CE608809668E24C61FC5BB8254BA36 -target:openmp -fusion:progress -scatter "e- e+ -> e-:e+ nue:nuebar u:ubar d:dbar"
!
! with all scattering amplitudes for the process(es)
!
!   flavor combinations:
!
!       1: e- e+ -> e+ nue ubar d
!       2: e- e+ -> e- nuebar u dbar
!
!   color flows:
!
!       1: (  0,  0) (  0,  0) -> (  0,  0) (  0,  0) (  0, -1) (  1,  0)
!       2: (  0,  0) (  0,  0) -> (  0,  0) (  0,  0) (  1,  0) (  0, -1)
!
!     NB: i.g. not all color flows contribute to all flavor
!     combinations.  Consult the array FLV_COL_IS_ALLOWED
!     below for the allowed combinations.
!
!   Color Factors:
!
!     (  1,  1): + N
!     (  2,  2): + N
!
!   vanishing or redundant flavor combinations:
!
!
!
module opr_ww_i1
  use kinds
  use omega95
  use omega_color, OCF => omega_color_factor
  use parameters_SM_ac
  implicit none
  private
  public :: number_particles_in, number_particles_out, number_color_indices, &
    reset_helicity_selection, new_event, is_allowed, get_amplitude, &
    color_sum, external_masses, openmp_supported, number_spin_states, &
    spin_states, number_flavor_states, flavor_states, number_color_flows, &
    color_flows, number_color_factors, color_factors, init, final, &
    update_alpha_s, md5sum

  ! DON'T EVEN THINK of removing the following!
  ! If the compiler complains about undeclared
  ! or undefined variables, you are compiling
  ! against an incompatible omega95 module!
  integer, dimension(7), parameter, private :: require = &
    (/ omega_spinors_2010_01_A, omega_spinor_cpls_2010_01_A, &
       omega_vectors_2010_01_A, omega_polarizations_2010_01_A, &
       omega_couplings_2010_01_A, omega_color_2010_01_A, &
       omega_utils_2010_01_A /)

  integer, parameter :: n_prt = 6
  integer, parameter :: n_in = 2
  integer, parameter :: n_out = 4
  integer, parameter :: n_cflow = 2
  integer, parameter :: n_cindex = 2
  integer, parameter :: n_flv = 2
  integer, parameter :: n_hel = 64
  integer, parameter :: n_co = 0
  integer, parameter :: n_cop = 0

  ! NB: you MUST NOT change the value of N_ here!!!
  !     It is defined here for convenience only and must be
  !     compatible with hardcoded values in the amplitude!
  real(kind=default), parameter :: N_ = 3
  logical, parameter :: F = .false.
  logical, parameter :: T = .true.

  integer, dimension(n_co,n_cop), save, protected :: table_coupling_orders

  integer, dimension(n_prt,n_hel), save, protected :: table_spin_states
  data table_spin_states(:,   1) / -1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,   2) / -1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,   3) / -1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,   4) / -1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,   5) / -1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,   6) / -1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,   7) / -1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,   8) / -1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,   9) / -1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  10) / -1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  11) / -1, -1,  1, -1,  1, -1 /
  data table_spin_states(:,  12) / -1, -1,  1, -1,  1,  1 /
  data table_spin_states(:,  13) / -1, -1,  1,  1, -1, -1 /
  data table_spin_states(:,  14) / -1, -1,  1,  1, -1,  1 /
  data table_spin_states(:,  15) / -1, -1,  1,  1,  1, -1 /
  data table_spin_states(:,  16) / -1, -1,  1,  1,  1,  1 /
  data table_spin_states(:,  17) / -1,  1, -1, -1, -1, -1 /
  data table_spin_states(:,  18) / -1,  1, -1, -1, -1,  1 /
  data table_spin_states(:,  19) / -1,  1, -1, -1,  1, -1 /
  data table_spin_states(:,  20) / -1,  1, -1, -1,  1,  1 /
  data table_spin_states(:,  21) / -1,  1, -1,  1, -1, -1 /
  data table_spin_states(:,  22) / -1,  1, -1,  1, -1,  1 /
  data table_spin_states(:,  23) / -1,  1, -1,  1,  1, -1 /
  data table_spin_states(:,  24) / -1,  1, -1,  1,  1,  1 /
  data table_spin_states(:,  25) / -1,  1,  1, -1, -1, -1 /
  data table_spin_states(:,  26) / -1,  1,  1, -1, -1,  1 /
  data table_spin_states(:,  27) / -1,  1,  1, -1,  1, -1 /
  data table_spin_states(:,  28) / -1,  1,  1, -1,  1,  1 /
  data table_spin_states(:,  29) / -1,  1,  1,  1, -1, -1 /
  data table_spin_states(:,  30) / -1,  1,  1,  1, -1,  1 /
  data table_spin_states(:,  31) / -1,  1,  1,  1,  1, -1 /
  data table_spin_states(:,  32) / -1,  1,  1,  1,  1,  1 /
  data table_spin_states(:,  33) /  1, -1, -1, -1, -1, -1 /
  data table_spin_states(:,  34) /  1, -1, -1, -1, -1,  1 /
  data table_spin_states(:,  35) /  1, -1, -1, -1,  1, -1 /
  data table_spin_states(:,  36) /  1, -1, -1, -1,  1,  1 /
  data table_spin_states(:,  37) /  1, -1, -1,  1, -1, -1 /
  data table_spin_states(:,  38) /  1, -1, -1,  1, -1,  1 /
  data table_spin_states(:,  39) /  1, -1, -1,  1,  1, -1 /
  data table_spin_states(:,  40) /  1, -1, -1,  1,  1,  1 /
  data table_spin_states(:,  41) /  1, -1,  1, -1, -1, -1 /
  data table_spin_states(:,  42) /  1, -1,  1, -1, -1,  1 /
  data table_spin_states(:,  43) /  1, -1,  1, -1,  1, -1 /
  data table_spin_states(:,  44) /  1, -1,  1, -1,  1,  1 /
  data table_spin_states(:,  45) /  1, -1,  1,  1, -1, -1 /
  data table_spin_states(:,  46) /  1, -1,  1,  1, -1,  1 /
  data table_spin_states(:,  47) /  1, -1,  1,  1,  1, -1 /
  data table_spin_states(:,  48) /  1, -1,  1,  1,  1,  1 /
  data table_spin_states(:,  49) /  1,  1, -1, -1, -1, -1 /
  data table_spin_states(:,  50) /  1,  1, -1, -1, -1,  1 /
  data table_spin_states(:,  51) /  1,  1, -1, -1,  1, -1 /
  data table_spin_states(:,  52) /  1,  1, -1, -1,  1,  1 /
  data table_spin_states(:,  53) /  1,  1, -1,  1, -1, -1 /
  data table_spin_states(:,  54) /  1,  1, -1,  1, -1,  1 /
  data table_spin_states(:,  55) /  1,  1, -1,  1,  1, -1 /
  data table_spin_states(:,  56) /  1,  1, -1,  1,  1,  1 /
  data table_spin_states(:,  57) /  1,  1,  1, -1, -1, -1 /
  data table_spin_states(:,  58) /  1,  1,  1, -1, -1,  1 /
  data table_spin_states(:,  59) /  1,  1,  1, -1,  1, -1 /
  data table_spin_states(:,  60) /  1,  1,  1, -1,  1,  1 /
  data table_spin_states(:,  61) /  1,  1,  1,  1, -1, -1 /
  data table_spin_states(:,  62) /  1,  1,  1,  1, -1,  1 /
  data table_spin_states(:,  63) /  1,  1,  1,  1,  1, -1 /
  data table_spin_states(:,  64) /  1,  1,  1,  1,  1,  1 /

  integer, dimension(n_prt,n_flv), save, protected :: table_flavor_states
  data table_flavor_states(:,   1) /  11, -11, -11,  12,  -2,   1 / ! e- e+ e+ nue ubar d
  data table_flavor_states(:,   2) /  11, -11,  11, -12,   2,  -1 / ! e- e+ e- nuebar u dbar

  integer, dimension(n_cindex,n_prt,n_cflow), save, protected :: table_color_flows
  data table_color_flows(:,:,   1) / 0,0,  0,0,  0,0,  0,0,  0,-1,  1,0 /
  data table_color_flows(:,:,   2) / 0,0,  0,0,  0,0,  0,0,  1,0,  0,-1 /

  logical, dimension(n_prt,n_cflow), save, protected :: table_ghost_flags
  data table_ghost_flags(:,   1) / F,  F,  F,  F,  F,  F /
  data table_ghost_flags(:,   2) / F,  F,  F,  F,  F,  F /

  integer, parameter :: n_cfactors = 2
  type(OCF), dimension(n_cfactors), save, protected :: table_color_factors
  real(kind=default), parameter, private :: color_factor_000001 = +N_
  data table_color_factors(     1) / OCF(1,1,color_factor_000001) /
  real(kind=default), parameter, private :: color_factor_000002 = +N_
  data table_color_factors(     2) / OCF(2,2,color_factor_000002) /

  logical, dimension(n_flv, n_cflow), save, protected ::  flv_col_is_allowed
  data flv_col_is_allowed(:,   1) / T, F /
  data flv_col_is_allowed(:,   2) / F, T /

  complex(kind=default), dimension(n_flv, n_cflow, n_hel), save :: amp

  logical, dimension(n_hel), save :: hel_is_allowed = T
  real(kind=default), dimension(n_hel), save :: hel_max_abs = 0
  real(kind=default), save :: hel_sum_abs = 0, hel_threshold = 1E10_default
  integer, save :: hel_count = 0, hel_cutoff = 100
  integer :: i
  integer, save, dimension(n_hel) :: hel_map = (/(i, i = 1, n_hel)/)
  integer, save :: hel_finite = n_hel

    type(momentum) :: p1, p2, p3, p4, p5, p6
    type(momentum) :: p12, p123, p124, p125, p126, p13, p134, p135, p136, &
      p14, p145, p146, p156, p23, p24, p34, p56
  type thread_local_data
    type(spinor) :: owf_fd1_i1_p6, owf_fu1_i1_p5, owf_fn1_p4, owf_fl1_p3, &
      owf_fl1_p1
    type(conjspinor) :: owf_fd1b_o1_p6, owf_fu1b_o1_p5, owf_fn1b_p4, &
      owf_fl1b_p3, owf_fl1b_p2
    type(spinor) :: owf_fd1_i1_p136, owf_fd1_i1_p126, owf_fd1_i1_p145, &
      owf_fu1_i1_p125, owf_fn1_p156, owf_fn1_p134, owf_fn1_p124, owf_fl1_p123
    type(conjspinor) :: owf_fd1b_o1_p126, owf_fu1b_o1_p146, owf_fu1b_o1_p135, &
      owf_fu1b_o1_p125, owf_fn1b_p124, owf_fl1b_p123
    type(vector) :: owf_fa_p23, owf_fa_p13, owf_fa_p12, owf_fz_p23, &
      owf_fz_p13, owf_fz_p12, owf_fwm_p56, owf_fwm_p34, owf_fwm_p14, &
      owf_fwp_p56, owf_fwp_p34, owf_fwp_p24
    complex(kind=default) :: oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1, &
      oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1
  end type thread_local_data
  type(thread_local_data) :: tld

contains

  pure function md5sum ()
    character(len=32) :: md5sum
    ! DON'T EVEN THINK of modifying the following line!
    md5sum = "84CE608809668E24C61FC5BB8254BA36"
  end function md5sum

  subroutine init (par, scheme)
    real(kind=default), dimension(*), intent(in) :: par
    integer, intent(in) :: scheme
    call import_from_whizard (par, scheme)
  end subroutine init

  subroutine final ()
  end subroutine final

  subroutine update_alpha_s (alpha_s)
    real(kind=default), intent(in) :: alpha_s
    call model_update_alpha_s (alpha_s)
  end subroutine update_alpha_s

  pure function number_particles_in () result (n)
    integer :: n
    n = n_in
  end function number_particles_in

  pure function number_particles_out () result (n)
    integer :: n
    n = n_out
  end function number_particles_out

  pure function number_spin_states () result (n)
    integer :: n
    n = size (table_spin_states, dim=2)
  end function number_spin_states

  pure subroutine spin_states (a)
    integer, dimension(:,:), intent(out) :: a
    a = table_spin_states
  end subroutine spin_states

  pure function number_flavor_states () result (n)
    integer :: n
    n = size (table_flavor_states, dim=2)
  end function number_flavor_states

  pure subroutine flavor_states (a)
    integer, dimension(:,:), intent(out) :: a
    a = table_flavor_states
  end subroutine flavor_states

  pure subroutine external_masses (m, flv)
    real(kind=default), dimension(:), intent(out) :: m
    integer, intent(in) :: flv
    select case (flv)
    case (  1)
      m( 1) = mass(11)
      m( 2) = mass(11)
      m( 3) = mass(11)
      m( 4) = mass(12)
      m( 5) = mass(2)
      m( 6) = mass(1)
    case (  2)
      m( 1) = mass(11)
      m( 2) = mass(11)
      m( 3) = mass(11)
      m( 4) = mass(12)
      m( 5) = mass(2)
      m( 6) = mass(1)
    end select
  end subroutine external_masses

  pure function openmp_supported () result (status)
    logical :: status
    status = .true.
  end function openmp_supported

  pure function number_color_indices () result (n)
    integer :: n
    n = size (table_color_flows, dim=1)
  end function number_color_indices

  pure function number_color_flows () result (n)
    integer :: n
    n = size (table_color_flows, dim=3)
  end function number_color_flows

  pure subroutine color_flows (a, g)
    integer, dimension(:,:,:), intent(out) :: a
    logical, dimension(:,:), intent(out) :: g
    a = table_color_flows
    g = table_ghost_flags
  end subroutine color_flows

  pure function number_color_factors () result (n)
    integer :: n
    n = size (table_color_factors)
  end function number_color_factors

  pure subroutine color_factors (cf)
    type(OCF), dimension(:), intent(out) :: cf
    cf = table_color_factors
  end subroutine color_factors

  function color_sum (flv, hel) result (amp2)
    integer, intent(in) :: flv, hel
    real(kind=default) :: amp2
    amp2 = real (omega_color_sum (flv, hel, amp, table_color_factors))
  end function color_sum

  subroutine new_event (p)
    real(kind=default), dimension(0:3,*), intent(in) :: p
    logical :: mask_dirty
    integer :: hel
    call calculate_amplitudes (amp, p, hel_is_allowed)
    if ((hel_threshold .gt. 0) .and. (hel_count .le. hel_cutoff)) then
      call omega_update_helicity_selection (hel_count, amp, hel_max_abs, &
              hel_sum_abs, hel_is_allowed, hel_threshold, hel_cutoff, &
              mask_dirty)
      if (mask_dirty) then
        hel_finite = 0
        do hel = 1, n_hel
          if (hel_is_allowed(hel)) then
            hel_finite = hel_finite + 1
            hel_map(hel_finite) = hel
          end if
        end do
      end if
    end if
  end subroutine new_event

  subroutine reset_helicity_selection (threshold, cutoff)
    real(kind=default), intent(in) :: threshold
    integer, intent(in) :: cutoff
    integer :: i
    hel_is_allowed = T
    hel_max_abs = 0
    hel_sum_abs = 0
    hel_count = 0
    hel_threshold = threshold
    hel_cutoff = cutoff
    hel_map = (/(i, i = 1, n_hel)/)
    hel_finite = n_hel
  end subroutine reset_helicity_selection

  pure function is_allowed (flv, hel, col) result (yorn)
    logical :: yorn
    integer, intent(in) :: flv, hel, col
    yorn = hel_is_allowed(hel) .and. flv_col_is_allowed(flv,col)
  end function is_allowed

  pure function get_amplitude (flv, hel, col) result (amp_result)
    complex(kind=default) :: amp_result
    integer, intent(in) :: flv, hel, col
    amp_result = amp(flv, col, hel)
  end function get_amplitude



  subroutine calculate_amplitudes (amp, k, mask)
    complex(kind=default), dimension(:,:,:), intent(out) :: amp
    real(kind=default), dimension(0:3,*), intent(in) :: k
    logical, dimension(:), intent(in) :: mask
    integer, dimension(n_prt) :: s
    integer :: h, hi
    p1 = - k(:,1) ! incoming
    p2 = - k(:,2) ! incoming
    p3 =   k(:,3) ! outgoing
    p4 =   k(:,4) ! outgoing
    p5 =   k(:,5) ! outgoing
    p6 =   k(:,6) ! outgoing
    p12 = p1 + p2
    p23 = p2 + p3
    p14 = p1 + p4
    p34 = p3 + p4
    p56 = p5 + p6
    p123 = p1 + p23
    p124 = p2 + p14
    p125 = p5 + p12
    p145 = p5 + p14
    p126 = p6 + p12
    p146 = p6 + p14
    p156 = p1 + p56
    p13 = p1 + p3
    p24 = p2 + p4
    p134 = p1 + p34
    p135 = p5 + p13
    p136 = p6 + p13
    amp = 0
    if (hel_finite == 0) return
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(s, h, tld) SCHEDULE(STATIC)
    do hi = 1, hel_finite
      h = hel_map(hi)
      s = table_spin_states(:,h)
      tld%owf_fl1_p1 = u (mass(11), - p1, s(1))
      tld%owf_fl1b_p2 = vbar (mass(11), - p2, s(2))
      tld%owf_fl1_p3 = v (mass(11), p3, s(3))
      tld%owf_fn1b_p4 = ubar (mass(12), p4, s(4))
      tld%owf_fu1_i1_p5 = v (mass(2), p5, s(5))
      tld%owf_fd1b_o1_p6 = ubar (mass(1), p6, s(6))
      tld%owf_fl1b_p3 = ubar (mass(11), p3, s(3))
      tld%owf_fn1_p4 = v (mass(12), p4, s(4))
      tld%owf_fu1b_o1_p5 = ubar (mass(2), p5, s(5))
      tld%owf_fd1_i1_p6 = v (mass(1), p6, s(6))
      call compute_fusions_0001 (tld)
      call compute_fusions_0002 (tld)
      call compute_fusions_0003 (tld)
      call compute_brakets_0001 (tld)
      call compute_brakets_0002 (tld)
      amp(1,1,h) = tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1
      amp(2,2,h) = tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1
    end do
!$OMP END PARALLEL DO
  end subroutine calculate_amplitudes
  subroutine compute_fusions_0001 (tld)
  type(thread_local_data), intent(inout) :: tld
      tld%owf_fa_p12 = pr_feynman(p12, &
         + v_ff(qlep,tld%owf_fl1b_p2,tld%owf_fl1_p1))
      tld%owf_fz_p12 = &
         pr_unitarity(p12,mass(23),wd_tl(p12,width(23)),.false., &
         + va_ff(gnclep(1),gnclep(2),tld%owf_fl1b_p2,tld%owf_fl1_p1))
      tld%owf_fa_p23 = pr_feynman(p23, &
         + v_ff(qlep,tld%owf_fl1b_p2,tld%owf_fl1_p3))
      tld%owf_fz_p23 = &
         pr_unitarity(p23,mass(23),wd_tl(p23,width(23)),.false., &
         + va_ff(gnclep(1),gnclep(2),tld%owf_fl1b_p2,tld%owf_fl1_p3))
      tld%owf_fwm_p14 = &
         pr_unitarity(p14,mass(24),wd_tl(p14,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fn1b_p4,tld%owf_fl1_p1))
      tld%owf_fwm_p34 = &
         pr_unitarity(p34,mass(24),wd_tl(p34,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fn1b_p4,tld%owf_fl1_p3))
      tld%owf_fwp_p56 = &
         pr_unitarity(p56,mass(24),wd_tl(p56,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fd1b_o1_p6,tld%owf_fu1_i1_p5))
      tld%owf_fl1_p123 = pr_psi(p123,mass(11),wd_tl(p123,width(11)),.false., &
         + f_vf(qlep,tld%owf_fa_p23,tld%owf_fl1_p1) &
         + f_vaf(gnclep(1),gnclep(2),tld%owf_fz_p23,tld%owf_fl1_p1) &
         - f_vf(qlep,tld%owf_fa_p12,tld%owf_fl1_p3) &
         - f_vaf(gnclep(1),gnclep(2),tld%owf_fz_p12,tld%owf_fl1_p3))
      tld%owf_fn1b_p124 = &
         pr_psibar(p124,mass(12),wd_tl(p124,width(12)),.false., &
         + f_fvl(gcc,tld%owf_fl1b_p2,tld%owf_fwm_p14) &
         - f_fva(gncneu(1),gncneu(2),tld%owf_fn1b_p4,tld%owf_fz_p12))
      tld%owf_fu1_i1_p125 = pr_psi(p125,mass(2),wd_tl(p125,width(2)),.false., &
         - f_vf(qup,tld%owf_fa_p12,tld%owf_fu1_i1_p5) &
         - f_vaf(gncup(1),gncup(2),tld%owf_fz_p12,tld%owf_fu1_i1_p5))
  end subroutine compute_fusions_0001
  subroutine compute_fusions_0002 (tld)
  type(thread_local_data), intent(inout) :: tld
      tld%owf_fd1_i1_p145 = pr_psi(p145,mass(1),wd_tl(p145,width(1)),.false., &
         - f_vlf(gcc,tld%owf_fwm_p14,tld%owf_fu1_i1_p5))
      tld%owf_fd1b_o1_p126 = &
         pr_psibar(p126,mass(1),wd_tl(p126,width(1)),.false., &
         - f_fv(qdwn,tld%owf_fd1b_o1_p6,tld%owf_fa_p12) &
         - f_fva(gncdwn(1),gncdwn(2),tld%owf_fd1b_o1_p6,tld%owf_fz_p12))
      tld%owf_fu1b_o1_p146 = &
         pr_psibar(p146,mass(2),wd_tl(p146,width(2)),.false., &
         - f_fvl(gcc,tld%owf_fd1b_o1_p6,tld%owf_fwm_p14))
      tld%owf_fn1_p156 = pr_psi(p156,mass(12),wd_tl(p156,width(12)),.false., &
         + f_vlf(gcc,tld%owf_fwp_p56,tld%owf_fl1_p1))
      tld%owf_fa_p13 = pr_feynman(p13, &
         + v_ff(qlep,tld%owf_fl1b_p3,tld%owf_fl1_p1))
      tld%owf_fz_p13 = &
         pr_unitarity(p13,mass(23),wd_tl(p13,width(23)),.false., &
         + va_ff(gnclep(1),gnclep(2),tld%owf_fl1b_p3,tld%owf_fl1_p1))
      tld%owf_fwp_p24 = &
         pr_unitarity(p24,mass(24),wd_tl(p24,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fl1b_p2,tld%owf_fn1_p4))
      tld%owf_fwp_p34 = &
         pr_unitarity(p34,mass(24),wd_tl(p34,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fl1b_p3,tld%owf_fn1_p4))
      tld%owf_fwm_p56 = &
         pr_unitarity(p56,mass(24),wd_tl(p56,width(24)),.false., &
         + vl_ff(gcc,tld%owf_fu1b_o1_p5,tld%owf_fd1_i1_p6))
      tld%owf_fl1b_p123 = &
         pr_psibar(p123,mass(11),wd_tl(p123,width(11)),.false., &
         + f_fv(qlep,tld%owf_fl1b_p2,tld%owf_fa_p13) &
         + f_fva(gnclep(1),gnclep(2),tld%owf_fl1b_p2,tld%owf_fz_p13) &
         - f_fv(qlep,tld%owf_fl1b_p3,tld%owf_fa_p12) &
         - f_fva(gnclep(1),gnclep(2),tld%owf_fl1b_p3,tld%owf_fz_p12))
  end subroutine compute_fusions_0002
  subroutine compute_fusions_0003 (tld)
  type(thread_local_data), intent(inout) :: tld
      tld%owf_fn1_p124 = pr_psi(p124,mass(12),wd_tl(p124,width(12)),.false., &
         + f_vlf(gcc,tld%owf_fwp_p24,tld%owf_fl1_p1) &
         - f_vaf(gncneu(1),gncneu(2),tld%owf_fz_p12,tld%owf_fn1_p4))
      tld%owf_fn1_p134 = pr_psi(p134,mass(12),wd_tl(p134,width(12)),.false., &
         + f_vlf(gcc,tld%owf_fwp_p34,tld%owf_fl1_p1) &
         - f_vaf(gncneu(1),gncneu(2),tld%owf_fz_p13,tld%owf_fn1_p4))
      tld%owf_fu1b_o1_p125 = &
         pr_psibar(p125,mass(2),wd_tl(p125,width(2)),.false., &
         - f_fv(qup,tld%owf_fu1b_o1_p5,tld%owf_fa_p12) &
         - f_fva(gncup(1),gncup(2),tld%owf_fu1b_o1_p5,tld%owf_fz_p12))
      tld%owf_fu1b_o1_p135 = &
         pr_psibar(p135,mass(2),wd_tl(p135,width(2)),.false., &
         - f_fv(qup,tld%owf_fu1b_o1_p5,tld%owf_fa_p13) &
         - f_fva(gncup(1),gncup(2),tld%owf_fu1b_o1_p5,tld%owf_fz_p13))
      tld%owf_fd1_i1_p126 = pr_psi(p126,mass(1),wd_tl(p126,width(1)),.false., &
         - f_vf(qdwn,tld%owf_fa_p12,tld%owf_fd1_i1_p6) &
         - f_vaf(gncdwn(1),gncdwn(2),tld%owf_fz_p12,tld%owf_fd1_i1_p6))
      tld%owf_fd1_i1_p136 = pr_psi(p136,mass(1),wd_tl(p136,width(1)),.false., &
         - f_vf(qdwn,tld%owf_fa_p13,tld%owf_fd1_i1_p6) &
         - f_vaf(gncdwn(1),gncdwn(2),tld%owf_fz_p13,tld%owf_fd1_i1_p6))
  end subroutine compute_fusions_0003
  subroutine compute_brakets_0001 (tld)
  type(thread_local_data), intent(inout) :: tld
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = 0
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + ( &
         - f_fvl(gcc,tld%owf_fn1b_p4,tld%owf_fwp_p56))*tld%owf_fl1_p123
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fn1b_p124*( &
         - f_vlf(gcc,tld%owf_fwp_p56,tld%owf_fl1_p3))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + ( &
         - f_fva(gncneu(1),gncneu(2),tld%owf_fn1b_p4,tld%owf_fz_p23) &
         + f_fvl(gcc,tld%owf_fl1b_p2,tld%owf_fwm_p34))*tld%owf_fn1_p156
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fu1b_o1_p146*( &
         + f_vaf(gncup(1),gncup(2),tld%owf_fz_p23,tld%owf_fu1_i1_p5) &
         + f_vf(qup,tld%owf_fa_p23,tld%owf_fu1_i1_p5))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + ( &
         - f_fvl(gcc,tld%owf_fd1b_o1_p6,tld%owf_fwm_p34))*tld%owf_fu1_i1_p125
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fd1b_o1_p126*( &
         - f_vlf(gcc,tld%owf_fwm_p34,tld%owf_fu1_i1_p5))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + ( &
         + f_fva(gncdwn(1),gncdwn(2),tld%owf_fd1b_o1_p6,tld%owf_fz_p23) &
         + f_fv(qdwn,tld%owf_fd1b_o1_p6,tld%owf_fa_p23))*tld%owf_fd1_i1_p145
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fa_p12*( &
         + tkv_vv((-ig1a),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + (-1)*tv_kvv((-ig1pkpg4a),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34) &
         + (-1)*tv_kvv(ig1pkmg4a,tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + lv_kvv((-ig1mkpg4a),tld%owf_fwm_p34,p34,tld%owf_fwp_p56) &
         + lv_kvv(ig1mkmg4a,tld%owf_fwp_p56,p56,tld%owf_fwm_p34) &
         + t5kv_vv((-rg5a),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + l5kv_vv((-ik5a),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + kg_kgkg((-ila),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34) &
         + kg5_kgkg((-il5a),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fwm_p14*( &
         - (-1)*tv_kvv((-ig1a),tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - (-1)*tv_kvv((-ig1z),tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56) &
         - tv_kvv((-ig1pkpg4a),tld%owf_fwp_p56,p56,tld%owf_fa_p23,p23) &
         - tv_kvv((-ig1pkpg4z),tld%owf_fwp_p56,p56,tld%owf_fz_p23,p23) &
         - tkv_vv(ig1pkmg4a,tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - tkv_vv(ig1pkmg4z,tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56) &
         - lkv_vv((-ig1mkpg4a),tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - lkv_vv((-ig1mkpg4z),tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56) &
         - lv_kvv(ig1mkmg4a,tld%owf_fwp_p56,p56,tld%owf_fa_p23) &
         - lv_kvv(ig1mkmg4z,tld%owf_fwp_p56,p56,tld%owf_fz_p23) &
         - t5v_kvv((-rg5a),tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - t5v_kvv((-rg5z),tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56) &
         - (-1)*l5v_kvv((-ik5a),tld%owf_fa_p23,p23,tld%owf_fwp_p56) &
         - (-1)*l5v_kvv((-ik5z),tld%owf_fz_p23,p23,tld%owf_fwp_p56) &
         - kg_kgkg((-ila),tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - kg_kgkg((-ilz),tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56) &
         - kg_kg5kg((-il5a),tld%owf_fa_p23,p23,tld%owf_fwp_p56,p56) &
         - kg_kg5kg((-il5z),tld%owf_fz_p23,p23,tld%owf_fwp_p56,p56))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
        tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 + tld%owf_fz_p12*( &
         + tkv_vv((-ig1z),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + (-1)*tv_kvv((-ig1pkpg4z),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34) &
         + (-1)*tv_kvv(ig1pkmg4z,tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + lv_kvv((-ig1mkpg4z),tld%owf_fwm_p34,p34,tld%owf_fwp_p56) &
         + lv_kvv(ig1mkmg4z,tld%owf_fwp_p56,p56,tld%owf_fwm_p34) &
         + t5kv_vv((-rg5z),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + l5kv_vv((-ik5z),tld%owf_fwm_p34,p34,tld%owf_fwp_p56,p56) &
         + kg_kgkg((-ilz),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34) &
         + kg5_kgkg((-il5z),tld%owf_fwp_p56,p56,tld%owf_fwm_p34,p34))
      tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 = &
         - tld%oks_fl1_fl1b_fl1b_fn1_fu1b_o1_fd1_i1 ! 4 vertices, 3 propagators
      ! unit symmetry factor
  end subroutine compute_brakets_0001
  subroutine compute_brakets_0002 (tld)
  type(thread_local_data), intent(inout) :: tld
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = 0
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fl1b_p123*( &
         - f_vlf(gcc,tld%owf_fwm_p56,tld%owf_fn1_p4))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + ( &
         + f_fvl(gcc,tld%owf_fl1b_p2,tld%owf_fwm_p56))*tld%owf_fn1_p134
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + ( &
         - f_fvl(gcc,tld%owf_fl1b_p3,tld%owf_fwm_p56))*tld%owf_fn1_p124
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fu1b_o1_p135*( &
         + f_vlf(gcc,tld%owf_fwp_p24,tld%owf_fd1_i1_p6))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fu1b_o1_p125*( &
         - f_vlf(gcc,tld%owf_fwp_p34,tld%owf_fd1_i1_p6))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + ( &
         + f_fvl(gcc,tld%owf_fu1b_o1_p5,tld%owf_fwp_p24))*tld%owf_fd1_i1_p136
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + ( &
         - f_fvl(gcc,tld%owf_fu1b_o1_p5,tld%owf_fwp_p34))*tld%owf_fd1_i1_p126
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fa_p12*( &
         + tkv_vv((-ig1a),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + (-1)*tv_kvv((-ig1pkpg4a),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56) &
         + (-1)*tv_kvv(ig1pkmg4a,tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + lv_kvv((-ig1mkpg4a),tld%owf_fwm_p56,p56,tld%owf_fwp_p34) &
         + lv_kvv(ig1mkmg4a,tld%owf_fwp_p34,p34,tld%owf_fwm_p56) &
         + t5kv_vv((-rg5a),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + l5kv_vv((-ik5a),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + kg_kgkg((-ila),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56) &
         + kg5_kgkg((-il5a),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fa_p13*( &
         - tkv_vv((-ig1a),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - (-1)*tv_kvv((-ig1pkpg4a),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56) &
         - (-1)*tv_kvv(ig1pkmg4a,tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - lv_kvv((-ig1mkpg4a),tld%owf_fwm_p56,p56,tld%owf_fwp_p24) &
         - lv_kvv(ig1mkmg4a,tld%owf_fwp_p24,p24,tld%owf_fwm_p56) &
         - t5kv_vv((-rg5a),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - l5kv_vv((-ik5a),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - kg_kgkg((-ila),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56) &
         - kg5_kgkg((-il5a),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fz_p12*( &
         + tkv_vv((-ig1z),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + (-1)*tv_kvv((-ig1pkpg4z),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56) &
         + (-1)*tv_kvv(ig1pkmg4z,tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + lv_kvv((-ig1mkpg4z),tld%owf_fwm_p56,p56,tld%owf_fwp_p34) &
         + lv_kvv(ig1mkmg4z,tld%owf_fwp_p34,p34,tld%owf_fwm_p56) &
         + t5kv_vv((-rg5z),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + l5kv_vv((-ik5z),tld%owf_fwm_p56,p56,tld%owf_fwp_p34,p34) &
         + kg_kgkg((-ilz),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56) &
         + kg5_kgkg((-il5z),tld%owf_fwp_p34,p34,tld%owf_fwm_p56,p56))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
        tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 + tld%owf_fz_p13*( &
         - tkv_vv((-ig1z),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - (-1)*tv_kvv((-ig1pkpg4z),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56) &
         - (-1)*tv_kvv(ig1pkmg4z,tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - lv_kvv((-ig1mkpg4z),tld%owf_fwm_p56,p56,tld%owf_fwp_p24) &
         - lv_kvv(ig1mkmg4z,tld%owf_fwp_p24,p24,tld%owf_fwm_p56) &
         - t5kv_vv((-rg5z),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - l5kv_vv((-ik5z),tld%owf_fwm_p56,p56,tld%owf_fwp_p24,p24) &
         - kg_kgkg((-ilz),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56) &
         - kg5_kgkg((-il5z),tld%owf_fwp_p24,p24,tld%owf_fwm_p56,p56))
      tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 = &
         - tld%oks_fl1_fl1b_fl1_fn1b_fu1_i1_fd1b_o1 ! 4 vertices, 3 propagators
      ! unit symmetry factor
  end subroutine compute_brakets_0002

end module opr_ww_i1
