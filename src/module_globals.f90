module module_globals
use,intrinsic :: iso_fortran_env,only:real32,real64
implicit none

  real(kind=real64), parameter :: dFillValue = -(HUGE(1.d0))
  real(kind=real32), parameter :: fFillValue = 1E-35
  real(kind=real32), parameter :: missing = fFillValue
  integer          , parameter :: iFillValue = -(HUGE(1))
  real(kind=real64), parameter :: g = 9.80665D0, rg = 1.D0/g
  real(kind=real32), parameter :: r_d = 287      ! [J/kg/K] gas const. dry air
  real(kind=real32), parameter :: cp = 1004.5    ! [J/kg/K] specific heat of dry air 
  real(kind=real32), parameter :: p1000mb = 100000
  real(kind=real32), parameter :: pi = 3.14159
  character(len=80), parameter :: namelist_name = "namelist.csi"

  character(len=120)  :: file_in, file_out
  integer  :: ncid_in, ncid_out
  ! integer  :: dims_w(4), dims_u(4), dims_v(4), dims_m(4)
  integer  :: ntimes

  real(kind=real32),dimension(:)  ,allocatable :: lon1d, lat1d, levels
  real(kind=real32),dimension(:,:),allocatable :: lon2d, lat2d


  integer :: v_id, u_id, w_id, p_id, ph_id 
  integer :: var_dimIDs(3) 
  integer :: xDimIDout=-1, yDimIDout=-1, zDimIDout=-1, tDimIDout=-1,    &
             xDimIDout_stag=-1, yDimIDout_stag=-1, zDimIDout_stag=-1
  integer :: xDimIDin=-1, yDimIDin=-1, zDimIDin=-1, tDimIDin=-1,        &
             xDimIDin_stag=-1, yDimIDin_stag=-1, zDimIDin_stag=-1
  integer :: xDimSIZE=-1, yDimSIZE=-1, zDimSIZE=-1, tDimSIZE=-1,        &
             xDimSIZE_stag=-1, yDimSIZE_stag=-1, zDimSIZE_stag=-1

  integer           :: nb_ticks_sec, nb_ticks_max, nb_ticks_initial
  real(kind=real32) :: elapsed_time

! NAMELIST SCHEME
  logical :: iswrf=.true.
  logical :: regional=.true.
  logical :: smooth=.false.
  real(kind=real32) :: smooth_sigma = 9.
  real(kind=real32) :: smooth_truncate = 3.
  logical :: output_height=.false.
  logical :: output_pressure=.false.
  logical :: debug=.false.
  logical :: debug_output=.false.
  logical :: q_criterion = .true.
  logical :: delta_criterion = .true.
  logical :: l2_criterion = .true.

  namelist /settings/ iswrf,regional,smooth,smooth_sigma,smooth_truncate, &
                      output_height,output_pressure,debug,debug_output
  namelist /criteria/ q_criterion,delta_criterion,l2_criterion

end module module_globals
