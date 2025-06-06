! WHIZARD matrix-element code interface
!
! Automatically generated file, do not edit

! Module: define library driver as an extension of the abstract driver type.
! This is used _only_ by the library dispatcher of a static executable.
! For a dynamical library, the stand-alone procedures are linked via libdl.

module default_lib_driver
  use iso_c_binding
  use iso_varying_string, string_t => varying_string
  use diagnostics
  use prclib_interfaces

  implicit none

  type, extends (prclib_driver_t) :: default_lib_driver_t
   contains
     procedure :: get_c_funptr => default_lib_driver_get_c_funptr
  end type default_lib_driver_t

contains

  function default_lib_driver_get_c_funptr (driver, feature) result (c_fptr)
    class(default_lib_driver_t), intent(inout) :: driver
    type(string_t), intent(in) :: feature
    type(c_funptr) :: c_fptr
    procedure(prc_get_n_processes) &
         :: default_lib_get_n_processes
    procedure(prc_get_stringptr) &
         :: default_lib_get_process_id_ptr
    procedure(prc_get_stringptr) &
         :: default_lib_get_model_name_ptr
    procedure(prc_get_stringptr) &
         :: default_lib_get_md5sum_ptr
    procedure(prc_get_log) &
         :: default_lib_get_openmp_status
    procedure(prc_get_int) &
         :: default_lib_get_n_in
    procedure(prc_get_int) &
         :: default_lib_get_n_out
    procedure(prc_get_int) &
         :: default_lib_get_n_flv
    procedure(prc_get_int) &
         :: default_lib_get_n_hel
    procedure(prc_get_int) &
         :: default_lib_get_n_col
    procedure(prc_get_int) &
         :: default_lib_get_n_cin
    procedure(prc_get_int) &
         :: default_lib_get_n_cf
    procedure(prc_set_int_tab1) &
         :: default_lib_set_flv_state_ptr
    procedure(prc_set_int_tab1) &
         :: default_lib_set_hel_state_ptr
    procedure(prc_set_col_state) &
         :: default_lib_set_col_state_ptr
    procedure(prc_set_color_factors) &
         :: default_lib_set_color_factors_ptr
    procedure(prc_get_fptr) &
         :: default_lib_get_fptr
    select case (char (feature))
    case ('get_n_processes')
       c_fptr = c_funloc (default_lib_get_n_processes)
    case ('get_process_id_ptr')
       c_fptr = c_funloc (default_lib_get_process_id_ptr)
    case ('get_model_name_ptr')
       c_fptr = c_funloc (default_lib_get_model_name_ptr)
    case ('get_md5sum_ptr')
       c_fptr = c_funloc (default_lib_get_md5sum_ptr)
    case ('get_openmp_status')
       c_fptr = c_funloc (default_lib_get_openmp_status)
    case ('get_n_in')
       c_fptr = c_funloc (default_lib_get_n_in)
    case ('get_n_out')
       c_fptr = c_funloc (default_lib_get_n_out)
    case ('get_n_flv')
       c_fptr = c_funloc (default_lib_get_n_flv)
    case ('get_n_hel')
       c_fptr = c_funloc (default_lib_get_n_hel)
    case ('get_n_col')
       c_fptr = c_funloc (default_lib_get_n_col)
    case ('get_n_cin')
       c_fptr = c_funloc (default_lib_get_n_cin)
    case ('get_n_cf')
       c_fptr = c_funloc (default_lib_get_n_cf)
    case ('set_flv_state_ptr')
       c_fptr = c_funloc (default_lib_set_flv_state_ptr)
    case ('set_hel_state_ptr')
       c_fptr = c_funloc (default_lib_set_hel_state_ptr)
    case ('set_col_state_ptr')
       c_fptr = c_funloc (default_lib_set_col_state_ptr)
    case ('set_color_factors_ptr')
       c_fptr = c_funloc (default_lib_set_color_factors_ptr)
    case ('get_fptr')
       c_fptr = c_funloc (default_lib_get_fptr)
    case default
       call msg_bug ('prclib2 driver setup: unknown function name')
    end select
  end function default_lib_driver_get_c_funptr

end module default_lib_driver

! Stand-alone external procedures: used for both static and dynamic linkage

! The MD5 sum of the library
function default_lib_md5sum () result (md5sum)
  implicit none
  character(32) :: md5sum
  md5sum = '906D7A57E3E57789066650303B7E31D6'
end function default_lib_md5sum

! Return the number of processes in this library
function default_lib_get_n_processes () result (n) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int) :: n
  n = 1
end function default_lib_get_n_processes

! Return the process ID of process #i (as a C pointer to a character array)
subroutine default_lib_get_process_id_ptr (i, cptr, len) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), intent(in) :: i
  type(c_ptr), intent(out) :: cptr
  integer(c_int), intent(out) :: len
  character(kind=c_char), dimension(:), allocatable, target, save :: a
  interface
     subroutine default_lib_string_to_array (string, a)
       use iso_c_binding
       implicit none
       character(*), intent(in) :: string
       character(kind=c_char), dimension(:), allocatable, intent(out) :: a
     end subroutine default_lib_string_to_array
  end interface
  select case (i)
  case (0);  if (allocated (a))  deallocate (a)
  case (1);  call default_lib_string_to_array ('ww_i1', a)
  end select
  if (allocated (a)) then
     cptr = c_loc (a)
     len = size (a)
  else
     cptr = c_null_ptr
     len = 0
  end if
end subroutine default_lib_get_process_id_ptr

! Return the model name for process #i (as a C pointer to a character array)
subroutine default_lib_get_model_name_ptr (i, cptr, len) bind(C)
  use iso_c_binding
  implicit none
  integer(c_int), intent(in) :: i
  type(c_ptr), intent(out) :: cptr
  integer(c_int), intent(out) :: len
  character(kind=c_char), dimension(:), allocatable, target, save :: a
  interface
     subroutine default_lib_string_to_array (string, a)
       use iso_c_binding
       implicit none
       character(*), intent(in) :: string
       character(kind=c_char), dimension(:), allocatable, intent(out) :: a
     end subroutine default_lib_string_to_array
  end interface
  select case (i)
  case (0);  if (allocated (a))  deallocate (a)
  case (1);  call default_lib_string_to_array ('SM_ac', a)
  end select
  if (allocated (a)) then
     cptr = c_loc (a)
     len = size (a)
  else
     cptr = c_null_ptr
     len = 0
  end if
end subroutine default_lib_get_model_name_ptr

! Return the MD5 sum for the process configuration (as a C pointer to a character array)
subroutine default_lib_get_md5sum_ptr (i, cptr, len) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_md5sum => md5sum
  implicit none
  interface
     function default_lib_md5sum () result (md5sum)
       character(32) :: md5sum
     end function default_lib_md5sum
  end interface
  integer(c_int), intent(in) :: i
  type(c_ptr), intent(out) :: cptr
  integer(c_int), intent(out) :: len
  character(kind=c_char), dimension(32), target, save :: md5sum
  select case (i)
  case (0)
     call copy (default_lib_md5sum ())
     cptr = c_loc (md5sum)
  case (1)
     call copy (ww_i1_md5sum ())
     cptr = c_loc (md5sum)
  case default
     cptr = c_null_ptr
  end select
  len = 32
contains
  subroutine copy (md5sum_tmp)
    character, dimension(32), intent(in) :: md5sum_tmp
    md5sum = md5sum_tmp
  end subroutine copy
end subroutine default_lib_get_md5sum_ptr

! Auxiliary: convert character string to array pointer
subroutine default_lib_string_to_array (string, a)
  use iso_c_binding
  implicit none
  character(*), intent(in) :: string
  character(kind=c_char), dimension(:), allocatable, intent(out) :: a
  integer :: i
  allocate (a (len (string)))
  do i = 1, size (a)
     a(i) = string(i:i)
  end do
end subroutine default_lib_string_to_array

! Return the OpenMP support status
function default_lib_get_openmp_status (i) result (openmp_status) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_openmp_supported => openmp_supported
  implicit none
  integer(c_int), intent(in) :: i
  logical(c_bool) :: openmp_status
  select case (i)
  case (1);  openmp_status = ww_i1_openmp_supported ()
  end select
end function default_lib_get_openmp_status

! Return the value of n_in
function default_lib_get_n_in (pid) result (n_in) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_in => number_particles_in
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_in
  select case (pid)
  case (1);  n_in = ww_i1_n_in ()
  end select
end function default_lib_get_n_in

! Return the value of n_out
function default_lib_get_n_out (pid) result (n_out) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_out => number_particles_out
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_out
  select case (pid)
  case (1);  n_out = ww_i1_n_out ()
  end select
end function default_lib_get_n_out

! Return the value of n_flv
function default_lib_get_n_flv (pid) result (n_flv) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_flv => number_flavor_states
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_flv
  select case (pid)
  case (1);  n_flv = ww_i1_n_flv ()
  end select
end function default_lib_get_n_flv

! Return the value of n_hel
function default_lib_get_n_hel (pid) result (n_hel) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_hel => number_spin_states
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_hel
  select case (pid)
  case (1);  n_hel = ww_i1_n_hel ()
  end select
end function default_lib_get_n_hel

! Return the value of n_col
function default_lib_get_n_col (pid) result (n_col) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_col => number_color_flows
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_col
  select case (pid)
  case (1);  n_col = ww_i1_n_col ()
  end select
end function default_lib_get_n_col

! Return the value of n_cin
function default_lib_get_n_cin (pid) result (n_cin) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_cin => number_color_indices
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_cin
  select case (pid)
  case (1);  n_cin = ww_i1_n_cin ()
  end select
end function default_lib_get_n_cin

! Return the value of n_cf
function default_lib_get_n_cf (pid) result (n_cf) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_n_cf => number_color_factors
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int) :: n_cf
  select case (pid)
  case (1);  n_cf = ww_i1_n_cf ()
  end select
end function default_lib_get_n_cf

! Set table: flv_state
subroutine default_lib_set_flv_state_ptr (pid, flv_state, shape) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_flv_state => flavor_states
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int), dimension(*), intent(out) :: flv_state
  integer(c_int), dimension(2), intent(in) :: shape
  integer, dimension(:,:), allocatable :: flv_state_tmp
  integer :: i, j
  select case (pid)
  case (1)
     allocate (flv_state_tmp (shape(1), shape(2)))
     call ww_i1_flv_state (flv_state_tmp)
     forall (i=1:shape(1), j=1:shape(2)) 
        flv_state(i + shape(1)*(j-1)) = flv_state_tmp(i,j)
     end forall
  end select
end subroutine default_lib_set_flv_state_ptr

! Set table: hel_state
subroutine default_lib_set_hel_state_ptr (pid, hel_state, shape) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_hel_state => spin_states
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int), dimension(*), intent(out) :: hel_state
  integer(c_int), dimension(2), intent(in) :: shape
  integer, dimension(:,:), allocatable :: hel_state_tmp
  integer :: i, j
  select case (pid)
  case (1)
     allocate (hel_state_tmp (shape(1), shape(2)))
     call ww_i1_hel_state (hel_state_tmp)
     forall (i=1:shape(1), j=1:shape(2)) 
        hel_state(i + shape(1)*(j-1)) = hel_state_tmp(i,j)
     end forall
  end select
end subroutine default_lib_set_hel_state_ptr

! Set tables: col_state, ghost_flag
subroutine default_lib_set_col_state_ptr (pid, col_state, ghost_flag, shape) bind(C)
  use iso_c_binding
  use opr_ww_i1, only: ww_i1_col_state => color_flows
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int), dimension(*), intent(out) :: col_state
  logical(c_bool), dimension(*), intent(out) :: ghost_flag
  integer(c_int), dimension(3), intent(in) :: shape
  integer, dimension(:,:,:), allocatable :: col_state_tmp
  logical, dimension(:,:), allocatable :: ghost_flag_tmp
  integer :: i, j, k
  select case (pid)
  case (1)
  allocate (col_state_tmp (shape(1), shape(2), shape(3)))
     allocate (ghost_flag_tmp (shape(2), shape(3)))
     call ww_i1_col_state (col_state_tmp, ghost_flag_tmp)
     forall (i = 1:shape(2), j = 1:shape(3))
        forall (k = 1:shape(1))
           col_state(k + shape(1) * (i + shape(2)*(j-1) - 1)) = col_state_tmp(k,i,j)
        end forall
        ghost_flag(i + shape(2)*(j-1)) = ghost_flag_tmp(i,j)
     end forall
  end select
end subroutine default_lib_set_col_state_ptr

! Set tables: color factors
subroutine default_lib_set_color_factors_ptr (pid, cf_index1, cf_index2, color_factors, shape) bind(C)
  use iso_c_binding
  use kinds
  use omega_color
  use opr_ww_i1, only: ww_i1_color_factors => color_factors
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int), dimension(1), intent(in) :: shape
  integer(c_int), dimension(*), intent(out) :: cf_index1, cf_index2
  complex(c_default_complex), dimension(*), intent(out) :: color_factors
  type(omega_color_factor), dimension(:), allocatable :: cf
  select case (pid)
  case (1)
     allocate (cf (shape(1)))
     call ww_i1_color_factors (cf)
     cf_index1(1:shape(1)) = cf%i1
     cf_index2(1:shape(1)) = cf%i2
     color_factors(1:shape(1)) = cf%factor
  end select
end subroutine default_lib_set_color_factors_ptr

! Return C pointer to a procedure:
! pid = process index;  fid = function index
subroutine default_lib_get_fptr (pid, fid, fptr) bind(C)
  use iso_c_binding
  use kinds
  implicit none
  integer(c_int), intent(in) :: pid
  integer(c_int), intent(in) :: fid
  type(c_funptr), intent(out) :: fptr
  interface
     subroutine ww_i1_init (par, scheme) bind(C)
       import
       real(c_default_float), dimension(*), intent(in) :: par
       integer(c_int), intent(in) :: scheme
     end subroutine ww_i1_init
  end interface
  interface
     subroutine ww_i1_update_alpha_s (alpha_s) bind(C)
       import
       real(c_default_float), intent(in) :: alpha_s
     end subroutine ww_i1_update_alpha_s
  end interface
  interface
     subroutine ww_i1_reset_helicity_selection (threshold, cutoff) bind(C)
       import
       real(c_default_float), intent(in) :: threshold
       integer(c_int), intent(in) :: cutoff
     end subroutine ww_i1_reset_helicity_selection
  end interface
  interface
     subroutine ww_i1_is_allowed (flv, hel, col, flag) bind(C)
       import
       integer(c_int), intent(in) :: flv, hel, col
       logical(c_bool), intent(out) :: flag
     end subroutine ww_i1_is_allowed
  end interface
  interface
     subroutine ww_i1_new_event (p) bind(C)
       import
       real(c_default_float), dimension(0:3,*), intent(in) :: p
     end subroutine ww_i1_new_event
  end interface
  interface
     subroutine ww_i1_get_amplitude (flv, hel, col, amp) bind(C)
       import
       integer(c_int), intent(in) :: flv, hel, col
       complex(c_default_complex), intent(out) :: amp
     end subroutine ww_i1_get_amplitude
  end interface
  select case (pid)
  case (1)
     select case (fid)
     case (1);  fptr = c_funloc (ww_i1_init)
     case (2);  fptr = c_funloc (ww_i1_update_alpha_s)
     case (3);  fptr = c_funloc (ww_i1_reset_helicity_selection)
     case (4);  fptr = c_funloc (ww_i1_is_allowed)
     case (5);  fptr = c_funloc (ww_i1_new_event)
     case (6);  fptr = c_funloc (ww_i1_get_amplitude)
     end select
  end select
end subroutine default_lib_get_fptr

subroutine ww_i1_init (par, scheme) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  real(c_default_float), dimension(*), intent(in) :: par
  integer(c_int), intent(in) :: scheme
  call init (par, scheme)
end subroutine ww_i1_init

subroutine ww_i1_update_alpha_s (alpha_s) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  real(c_default_float), intent(in) :: alpha_s
  call update_alpha_s (alpha_s)
end subroutine ww_i1_update_alpha_s

subroutine ww_i1_reset_helicity_selection (threshold, cutoff) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  real(c_default_float), intent(in) :: threshold
  integer(c_int), intent(in) :: cutoff
  call reset_helicity_selection (threshold, int (cutoff))
end subroutine ww_i1_reset_helicity_selection

subroutine ww_i1_is_allowed (flv, hel, col, flag) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  integer(c_int), intent(in) :: flv, hel, col
  logical(c_bool), intent(out) :: flag
  flag = is_allowed (int (flv), int (hel), int (col))
end subroutine ww_i1_is_allowed

subroutine ww_i1_new_event (p) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  real(c_default_float), dimension(0:3,*), intent(in) :: p
  call new_event (p)
end subroutine ww_i1_new_event

subroutine ww_i1_get_amplitude (flv, hel, col, amp) bind(C)
  use iso_c_binding
  use kinds
  use opr_ww_i1
  integer(c_int), intent(in) :: flv, hel, col
  complex(c_default_complex), intent(out) :: amp
  amp = get_amplitude (int (flv), int (hel), int (col))
end subroutine ww_i1_get_amplitude
