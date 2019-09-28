module module_io
!===============================================================================
!
! module_io: The module provides some procedures for 
!               input|output operations. Here we have everything 
!               related to NetCDF and the namelist.
!               
!===============================================================================

use,intrinsic :: iso_fortran_env,only:real32,real64
use netcdf
use module_globals
implicit none

private

public :: print_header 
public :: print_footer 
public :: read_namelist 
public :: open_src_nc 
public :: create_dst_nc 
public :: finalize 
public :: put_var_time 
public :: get_nc_var 
public :: wrf_user_unstagger 
public :: save_kernel 
public :: print_stat 

public :: time_init 
public :: time_start 
public :: time_end 

    interface put_var_time
        procedure put_var_xyt_time, put_var_xyzt_time, &
                                    put_var_xyzt_time_64
    end interface put_var_time

    interface get_nc_var
        procedure get_nc_var2d, get_nc_var3d, get_nc_var1d, get_nc_var0d
    end interface get_nc_var

    interface print_stat
        procedure print_stat_3d_32, print_stat_3d_64
    end interface print_stat

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine get_nc_var3d(ncid, varname, itime, output)
use netcdf
implicit none
character(len=*), intent(in)    :: varname
integer,intent(in)              :: itime,ncid
real(kind=real32),intent(inout) :: output(:,:,:)
integer                         :: VarId, xdim, ydim, zdim
real(kind=real64)               :: add_offset,scale_factor
integer(2),dimension(:,:,:),allocatable :: sarray

    xdim  = UBOUND(output,1)
    ydim  = UBOUND(output,2)
    zdim  = UBOUND(output,3)

    call check( nf90_inq_varid(ncid, trim(varname), VarId) )
    
    add_offset = get_offset(ncid,VarId)
    scale_factor = get_scale(ncid,VarId)
    
    if(add_offset .ne. 0) then
        allocate( sarray(xdim,ydim,zdim) )

        call check( nf90_get_var(ncid, VarId, sarray , start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )
        output = float(sarray)*scale_factor+add_offset

        deallocate( sarray )
    else
      call check( nf90_get_var(ncid, VarId, output, start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /)  ) )
    end if

end subroutine get_nc_var3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine get_nc_var2d(ncid, varname, itime, output)
use netcdf
implicit none
character(len=*), intent(in)    :: varname
integer,intent(in)              :: itime,ncid
real(kind=real32),intent(out)   :: output(:,:)
integer                         :: VarId, ydim, xdim
real(kind=real64)               :: add_offset,scale_factor
integer(2),dimension(:,:),allocatable :: sarray

    xdim  = UBOUND(output,1)
    ydim  = UBOUND(output,2)

    call check( nf90_inq_varid(ncid, trim(varname), VarId) )
    
    add_offset = get_offset(ncid,VarId)
    scale_factor = get_scale(ncid,VarId)
    
    if(add_offset .ne. 0) then
        allocate( sarray(xdim,ydim) )

        call check( nf90_get_var(ncid, VarId, sarray , start = (/ 1, 1, itime /), count = (/ xdim, ydim, 1 /)  ) )
        output = float(sarray)*scale_factor+add_offset

        deallocate( sarray )
    else
      call check( nf90_get_var(ncid, VarId, output, start = (/ 1, 1, itime /), count = (/ xdim, ydim, 1 /)  ) )
    end if

end subroutine get_nc_var2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine get_nc_var1d(ncid, varname, itime, output)
use netcdf
implicit none
character(len=*), intent(in)    :: varname
integer,intent(in)              :: itime,ncid
real(kind=real32),intent(out)   :: output(:)
integer                         :: VarId, ydim, xdim
real(kind=real64)               :: add_offset,scale_factor
integer(2),dimension(:),allocatable :: sarray

    xdim  = UBOUND(output,1)

    call check( nf90_inq_varid(ncid, trim(varname), VarId) )

    add_offset = get_offset(ncid,VarId)
    scale_factor = get_scale(ncid,VarId)
    
    if(add_offset .ne. 0) then
        allocate( sarray(xdim) )

        call check( nf90_get_var(ncid, VarId, sarray , start = (/ 1, itime /), count = (/ xdim, 1 /) ) )
        output = float(sarray)*scale_factor+add_offset

        deallocate( sarray )
    else
      call check( nf90_get_var(ncid, VarId, output, start = (/ 1, itime /), count = (/ xdim, 1 /)  ) )
    end if

end subroutine get_nc_var1d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine get_nc_var0d(ncid, varname, itime, output)
use netcdf
implicit none
character(len=*), intent(in)    :: varname
integer,intent(in)              :: itime,ncid
real(kind=real32),intent(out)   :: output
integer                         :: VarId, ydim, xdim
real(kind=real64)               :: add_offset,scale_factor
integer(2)                      :: sarray

    call check( nf90_inq_varid(ncid, trim(varname), VarId) )
    
    add_offset = get_offset(ncid,VarId)
    scale_factor = get_scale(ncid,VarId)

    ! call check( nf90_get_var(ncid, VarId, output, start = (/ itime /) ) )

    if(add_offset .ne. 0) then

        call check( nf90_get_var(ncid, VarId, sarray , start = (/ itime /) ) )
        output = float(sarray)*scale_factor+add_offset

    else
        call check( nf90_get_var(ncid, VarId, output, start = (/ itime /) ) )
    end if


end subroutine get_nc_var0d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine put_var_xyt_time(ncid,array,var_name,var_longname,var_units,itime,dimids)
implicit none
    integer,intent(in)              :: ncid, itime
    character(len=*), intent(in)    :: var_name,var_longname,var_units
    real(kind=real32)               :: array(:,:)
    integer                         :: VAR_ID, status, xdim, ydim, dimids(3)

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)

    status = nf90_inq_varid(ncid, var_name, VAR_ID) 

    if(status /= nf90_noerr) then ! variable do not exists
        call check(  nf90_redef(ncid) )
        call check( nf90_def_var(ncid, name = trim(var_name), xtype = NF90_FLOAT, dimids = dimids, varid = VAR_ID) )
        call check( nf90_put_att(ncid, VAR_ID, "long_name", trim(var_longname)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "units", trim(var_units)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "_FillValue", fFillValue) )
        call check( nf90_enddef(ncid_out) )
        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, itime /), count = (/ xdim, ydim, 1 /) ) )

    else ! variable does not exist => create

        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, itime /), count = (/ xdim, ydim, 1 /) ) )

    end if


end subroutine put_var_xyt_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine put_var_xyzt_time(ncid,array,var_name,var_longname,var_units,itime,dimids)
implicit none
    integer,intent(in)              :: ncid, itime
    character(len=*), intent(in)    :: var_name,var_longname,var_units
    real(kind=real32)               :: array(:,:,:)
    integer                         :: VAR_ID, status, xdim, ydim, zdim, dimids(4)

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)
    zdim  = UBOUND(array,3)

    status = nf90_inq_varid(ncid, var_name, VAR_ID) 

    if(status /= nf90_noerr) then ! variable do not exists
        call check(  nf90_redef(ncid) )
        call check( nf90_def_var(ncid, name = trim(var_name), xtype = NF90_FLOAT, dimids = dimids, varid = VAR_ID) )
        call check( nf90_put_att(ncid, VAR_ID, "long_name", trim(var_longname)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "units", trim(var_units)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "_FillValue", fFillValue) )
        call check( nf90_enddef(ncid_out) )
        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )

    else 

        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )

    end if


end subroutine put_var_xyzt_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine put_var_xyzt_time_64(ncid,array,var_name,var_longname,var_units,itime,dimids)
implicit none
    integer,intent(in)              :: ncid, itime
    character(len=*), intent(in)    :: var_name,var_longname,var_units
    real(kind=real64)               :: array(:,:,:)
    integer                         :: VAR_ID, status, xdim, ydim, zdim, dimids(4)

    xdim  = UBOUND(array,1)
    ydim  = UBOUND(array,2)
    zdim  = UBOUND(array,3)

    status = nf90_inq_varid(ncid, var_name, VAR_ID) 

    if(status /= nf90_noerr) then ! variable do not exists
        call check(  nf90_redef(ncid) )
        call check( nf90_def_var(ncid, name = trim(var_name), xtype = NF90_DOUBLE, dimids = dimids, varid = VAR_ID) )
        call check( nf90_put_att(ncid, VAR_ID, "long_name", trim(var_longname)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "units", trim(var_units)) ) 
        call check( nf90_put_att(ncid, VAR_ID, "_FillValue", dFillValue) )
        call check( nf90_enddef(ncid_out) )
        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )

    else 

        call check( nf90_put_var(ncid, VAR_ID, array, start = (/ 1, 1, 1, itime /), count = (/ xdim, ydim, zdim, 1 /) ) )

    end if


end subroutine put_var_xyzt_time_64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine wrf_user_unstagger(stag_array,mass_array)
! inperlolates staggered grid (Arakawa C) WRF on mass points
! simular to ncl function
! stag_array - staggered grid
! mass_array - output mass grid
implicit none
    integer     :: xdim,ydim,zdim,nxdim,nydim,nzdim
    integer     :: ii, jj, kk
    integer     :: ndims_stag, ndims_mass
    real(kind=real32), intent(in)    :: stag_array(:,:,:)
    real(kind=real32), intent(out)   :: mass_array(:,:,:)
    
    ndims_stag = SIZE(SHAPE(stag_array))
    ndims_mass = SIZE(SHAPE(mass_array))
    
    if(ndims_stag.NE.3 .OR. ndims_mass.NE.3) then
        print*, "wrf_user_unstagger: ERROR! array has to be 3D. STOP."
        call finalize
    end if
    
    xdim  = UBOUND(stag_array,1)
    ydim  = UBOUND(stag_array,2)
    zdim  = UBOUND(stag_array,3)
    nxdim = UBOUND(mass_array,1)
    nydim = UBOUND(mass_array,2)
    nzdim = UBOUND(mass_array,3)

    mass_array = fFillValue

    if(xdim.GT.nxdim) then
        ! print*,"x unstaggering!!!!!!!!!!!!"
        do jj = 1, nydim
        do kk = 1, nzdim
        do ii = 1, nxdim
            if( stag_array(ii,jj,kk).NE.fFillValue .AND. stag_array(ii+1,jj,kk).NE.fFillValue) then
                mass_array(ii,jj,kk) = 0.5 * ( stag_array(ii+1,jj,kk) + stag_array(ii,jj,kk)  )
            end if
        end do
        end do
        end do
        ! mass_array = 0.5 * ( stag_array(2:xdim,:,:) + stag_array(1:xdim-1,:,:)  )
    elseif(ydim.GT.nydim) then
        ! print*,"x unstaggering!!!!!!!!!!!!"
        do jj = 1, nydim
        do kk = 1, nzdim
        do ii = 1, nxdim
            if( stag_array(ii,jj,kk).NE.fFillValue .AND. stag_array(ii,jj+1,kk).NE.fFillValue) then
                mass_array(ii,jj,kk) = 0.5 * ( stag_array(ii,jj+1,kk) + stag_array(ii,jj,kk)  )
            end if
        end do
        end do
        end do
        ! mass_array = 0.5 * ( stag_array(:,2:ydim,:) + stag_array(:,1:ydim-1,:)  )
    elseif(zdim.GT.nzdim) then
        ! print*,"z unstaggering"
        do kk = 1, nzdim
        do ii = 1, nxdim
        do jj = 1, nydim
            if( stag_array(ii,jj,kk).NE.fFillValue .AND. stag_array(ii,jj,kk+1).NE.fFillValue) then
                mass_array(ii,jj,kk) = 0.5 * ( stag_array(ii,jj,kk+1) + stag_array(ii,jj,kk)  )
            end if
        end do
        end do
        end do
    else
        if( xdim.EQ.nxdim .AND. ydim.EQ.nydim .AND. zdim.EQ.nzdim) then
            print*,"wrf_user_unstagger: NOTE. No unstaggering required, as the input field is already on mass points."
            mass_array = stag_array
        end if
    
    end if
    
end subroutine wrf_user_unstagger


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine get_ndims(varname,dims)
implicit none
character(len=*), intent(in)    :: varname
integer,intent(out)             :: dims(:)
integer                         :: VarId, ndims, ii
integer, dimension(nf90_max_var_dims)   :: dimIDs   

    dims = 0
    VarId = -999
    
    call check( nf90_inq_varid(ncid_in, trim(varname), VarId) )
    call check( nf90_inquire_variable(ncid_in, varid=VarId, ndims = ndims, dimids = dimIDs) )
    do ii = 1, ndims
        call check( nf90_inquire_dimension(ncid_in, dimIDs(ii), len = dims(ii)) )
    end do

end subroutine get_ndims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine print_header
implicit none
  
  print*,"-> Files:"
  print*,"   Input file :      ",trim(file_in)
  print*,"   Output file: ",trim(file_out)

  print*,"-> Settings:"
  if(iswrf) then
    print*,"   Input data is WRF"
  else 
    print*,"   Input data is not WRF"
  end if

end subroutine print_header

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine time_init
use ifport
implicit none
  
    CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)

end subroutine time_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine time_start(ticks)
use ifport
implicit none
  integer,intent(out) :: ticks
  
    CALL SYSTEM_CLOCK(COUNT=ticks)

end subroutine time_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine time_end(ticks_start,sec)
use ifport
implicit none
  integer,intent(in) :: ticks_start
  real(kind=real32),intent(out) :: sec
  integer :: ticks_end,ticks
  
  CALL SYSTEM_CLOCK(COUNT=ticks_end)
  ticks = ticks_end - ticks_start
  IF (ticks_end < ticks_start) ticks = ticks + nb_ticks_max
  sec   = float(ticks) / nb_ticks_sec

end subroutine time_end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine print_footer
implicit none
  
  call time_end(nb_ticks_initial,elapsed_time)
  print'(a,a,a,f7.1)',"Done with ",trim(file_in), &
                      "; Elapsed time (min): ", elapsed_time/60.

end subroutine print_footer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine create_dst_nc
implicit none

  integer :: TimeId, unit_num

  ! delete previous file
  unit_num = 200 
  open(unit_num,file=trim(file_out),status='unknown')
  close(unit_num,status='delete')

  call check( nf90_create( path = file_out ,              &
              cmode = or(nf90_clobber,nf90_64bit_offset), &
              ncid=ncid_out) )
  
  call copy_all_dims(ncid_in,ncid_out)

  call copy_time(ncid_in,ncid_out)

  call copy_coordinates(ncid_in,ncid_out)

end subroutine create_dst_nc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine copy_all_dims(ncid_in, ncid_out)
implicit none
  integer,intent(in) :: ncid_in,ncid_out
  integer :: numDims,ii, status
  character(len = nf90_max_name) :: dim_name
  integer :: dim_len
  integer :: dim_id
  character(len = 20) :: tDim_names(2)=(/"Time","time"/)
  character(len = 20) :: xDim_names(3)=(/"west_east","lon","longitude"/)
  character(len = 20) :: yDim_names(3)=(/"south_north","lat","latitude"/)
  character(len = 20) :: zDim_names(2)=(/"bottom_top","level"/)
  character(len = 20) :: xDim_stag_names(1)=(/"west_east_stag"/)
  character(len = 20) :: yDim_stag_names(1)=(/"south_north_stag"/)
  character(len = 20) :: zDim_stag_names(1)=(/"bottom_top_stag"/)

  ! copy all dimensions
  call check( nf90_inquire(ncid_in, nDimensions=numDims))
  do ii = 1, numDims
    dim_id = ii
    call check( nf90_inquire_dimension(ncid_in, dimid=dim_id, name=dim_name, len=dim_len) )
    call check( nf90_def_dim(ncid_out, name=dim_name, len=dim_len, dimid=dim_id ) )
  end do

  call check( nf90_enddef(ncid_out) )


  ! Retrieve IDs of all dimensions in src and dst files
  ! SRC DIMENSIONS
  do ii = 1, size(tDim_names)
    status = nf90_inq_dimid(ncid_in, tDim_names(ii), tDimIDin)
    if(status .eq. nf90_noerr) then
      call check ( nf90_inquire_dimension(ncid_in, dimid=tDimIDin, len=tDimSIZE) )
      exit
    end if
  end do 

  do ii = 1, size(xDim_names)
    status = nf90_inq_dimid(ncid_in, xDim_names(ii), xDimIDin)
    if(status .eq. nf90_noerr) then
      call check ( nf90_inquire_dimension(ncid_in, dimid=xDimIDin, len=xDimSIZE) )
      exit
    end if
  end do 

  do ii = 1, size(yDim_names)
    status = nf90_inq_dimid(ncid_in, yDim_names(ii), yDimIDin)
    if(status .eq. nf90_noerr) then
      call check ( nf90_inquire_dimension(ncid_in, dimid=yDimIDin, len=yDimSIZE) )
      exit
    end if
  end do 

  do ii = 1, size(zDim_names)
    status = nf90_inq_dimid(ncid_in, zDim_names(ii), zDimIDin)
    if(status .eq. nf90_noerr) then
      call check ( nf90_inquire_dimension(ncid_in, dimid=zDimIDin, len=zDimSIZE) )
      exit
    end if
  end do 

  if (iswrf) then

    do ii = 1, size(xDim_stag_names)
      status = nf90_inq_dimid(ncid_in, xDim_stag_names(ii), xDimIDin_stag)
      if(status .eq. nf90_noerr) then
        call check ( nf90_inquire_dimension(ncid_in, dimid=xDimIDin_stag, len=xDimSIZE_stag) )
        exit
      end if
    end do 

    do ii = 1, size(yDim_stag_names)
      status = nf90_inq_dimid(ncid_in, yDim_stag_names(ii), yDimIDin_stag)
      if(status .eq. nf90_noerr) then
        call check ( nf90_inquire_dimension(ncid_in, dimid=yDimIDin_stag, len=yDimSIZE_stag) )
        exit
      end if
    end do 

    do ii = 1, size(zDim_stag_names)
      status = nf90_inq_dimid(ncid_in, zDim_stag_names(ii), zDimIDin_stag)
      if(status .eq. nf90_noerr) then
        call check ( nf90_inquire_dimension(ncid_in, dimid=zDimIDin_stag, len=zDimSIZE_stag) )
        exit
      end if
    end do 

  end if

  ! DST DIMENSIONS

  do ii = 1, size(tDim_names)
    status = nf90_inq_dimid(ncid_out, tDim_names(ii), tDimIDout)
    if(status .eq. nf90_noerr) exit
  end do 

  do ii = 1, size(xDim_names)
    status = nf90_inq_dimid(ncid_out, xDim_names(ii), xDimIDout)
    if(status .eq. nf90_noerr) exit
  end do 

  do ii = 1, size(yDim_names)
    status = nf90_inq_dimid(ncid_out, yDim_names(ii), yDimIDout)
    if(status .eq. nf90_noerr) exit
  end do 

  do ii = 1, size(zDim_names)
    status = nf90_inq_dimid(ncid_out, zDim_names(ii), zDimIDout)
    if(status .eq. nf90_noerr) exit
  end do 

  if (iswrf) then

    do ii = 1, size(xDim_stag_names)
      status = nf90_inq_dimid(ncid_out, xDim_stag_names(ii), xDimIDout_stag)
      if(status .eq. nf90_noerr) exit
    end do 

    do ii = 1, size(yDim_stag_names)
      status = nf90_inq_dimid(ncid_out, yDim_stag_names(ii), yDimIDout_stag)
      if(status .eq. nf90_noerr) exit
    end do 

    do ii = 1, size(zDim_stag_names)
      status = nf90_inq_dimid(ncid_out, zDim_stag_names(ii), zDimIDout_stag)
      if(status .eq. nf90_noerr) exit
    end do 

  end if

end subroutine copy_all_dims

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine copy_time(ncid_in, ncid_out)
implicit none
  integer,intent(in) :: ncid_in,ncid_out
  character(len = nf90_max_name) :: time_names(3)=(/"Times","XTIME","time"/)
  character(len = nf90_max_name) :: time_name
  character(len = nf90_max_name) :: nameAtt
  integer :: ii, ia, TimeInId, TimeOutId, xtype, numDims, numAtts, status
  integer, dimension(nf90_max_var_dims) :: dim_id_var
  
  character (len=19),dimension(tDimSIZE) :: time_str
  real(kind=real64),dimension(tDimSIZE)  :: time_dbl
  real(kind=real32),dimension(tDimSIZE)  :: time_flt
  integer,dimension(tDimSIZE)            :: time_int

  do ii = 1, size(time_names)

    time_name = time_names(ii)
    status = nf90_inq_varid(ncid_in, time_name, TimeInId)

    if(status .eq. nf90_noerr) then 

      
      call check( nf90_redef(ncid_out)  )

      call check( nf90_inquire_variable(ncid_in, TimeInId, ndims = numDims,   &
                    dimids = dim_id_var, nAtts = numAtts, xtype = xtype ) )
      call check( nf90_def_var(ncid_out, name = time_name, xtype = xtype,     &
                  dimids = (/ dim_id_var(:numDims) /), varid = TimeOutId ) )

      do ia = 1, numAtts
        call check ( nf90_inq_attname(ncid_in, TimeInId, ia, nameAtt) )
        call check ( nf90_copy_att(ncid_in, TimeInId, nameAtt, ncid_out, TimeOutId) )
      end do

      call check( nf90_enddef(ncid_out) )

      if(xtype.eq.2) then ! string time  
        call check( nf90_get_var(ncid_in,  TimeInId,  time_str ) )
        call check( nf90_put_var(ncid_out, TimeOutId, time_str ) )
      else if(xtype.eq.4) then
        call check( nf90_get_var(ncid_in,  TimeInId,  time_int ) )
        call check( nf90_put_var(ncid_out, TimeOutId, time_int ) )
      else if(xtype.eq.5) then
        call check( nf90_get_var(ncid_in,  TimeInId,  time_flt ) )
        call check( nf90_put_var(ncid_out, TimeOutId, time_flt ) )
      else if(xtype.eq.6) then
        call check( nf90_get_var(ncid_in,  TimeInId,  time_dbl ) )
        call check( nf90_put_var(ncid_out, TimeOutId, time_dbl ) )
      else
        call finalize
        stop "copy_time: UNKNOWN TIME FORMAT IN SRC FILE. STOP."
      end if

    
    end if

  end do

end subroutine copy_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine copy_coordinates(ncid_in, ncid_out)
implicit none
  integer,intent(in) :: ncid_in,ncid_out
  character(len = nf90_max_name) :: lat_names(3)=(/"XLAT","lat","latitude"/)
  character(len = nf90_max_name) :: lon_names(3)=(/"XLONG","lon","longitude"/)
  character(len = nf90_max_name) :: lev_names(1)=(/"level"/)

  character(len = nf90_max_name) :: var_name
  character(len = nf90_max_name) :: nameAtt
  integer :: ii, ia, VarInId, VarOutId, xtype, numDims, numAtts, status
  integer, dimension(nf90_max_var_dims) :: dim_id_var
  integer, dimension(nf90_max_var_dims) :: dims
  

  do ii = 1, size(lat_names)

    var_name = lat_names(ii)
    status = nf90_inq_varid(ncid_in, var_name, VarInId)

    if(status .eq. nf90_noerr) then 
      
      call check( nf90_redef(ncid_out)  )

      call check( nf90_inquire_variable(ncid_in, VarInId, ndims = numDims,   &
                    dimids = dim_id_var, nAtts = numAtts, xtype = xtype ) )
      ! in WRF coordinates are written in 3D (moving nests accounted)
      ! We do not use moving nests, so converting 3D->2D
      if(numDims.ne.3) then
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:numDims) /), varid = VarOutId ) )
      else
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:2) /), varid = VarOutId ) )
      end if

      do ia = 1, numAtts
        call check ( nf90_inq_attname(ncid_in, VarInId, ia, nameAtt) )
        call check ( nf90_copy_att(ncid_in, VarInId, nameAtt, ncid_out, VarOutId) )
      end do

      call check( nf90_enddef(ncid_out) )

      call get_ndims(var_name,dims)

      if(numDims.eq.1) then
        allocate( lat1d( dims(1) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lat1d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lat1d ) )
      else if(numDims.eq.2) then
        allocate( lat2d( dims(1),dims(2) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lat2d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lat2d ) )
      else if(numDims.eq.3) then
        allocate( lat2d( dims(1),dims(2) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lat2d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lat2d ) )
      else
        call finalize
        stop "copy_coordinates: UNKNOWN LAT FORMAT IN SRC FILE. STOP."
      end if

    
     end if

  end do

  ! Longitude
  do ii = 1, size(lon_names)

    var_name = lon_names(ii)
    status = nf90_inq_varid(ncid_in, var_name, VarInId)

    if(status .eq. nf90_noerr) then 
      
      call check( nf90_redef(ncid_out)  )

      call check( nf90_inquire_variable(ncid_in, VarInId, ndims = numDims,   &
                    dimids = dim_id_var, nAtts = numAtts, xtype = xtype ) )
      ! in WRF coordinates are written in 3D (moving nests accounted)
      ! We do not use moving nests, so converting 3D->2D
      if(numDims.ne.3) then
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:numDims) /), varid = VarOutId ) )
      else
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:2) /), varid = VarOutId ) )
      end if

      do ia = 1, numAtts
        call check ( nf90_inq_attname(ncid_in, VarInId, ia, nameAtt) )
        call check ( nf90_copy_att(ncid_in, VarInId, nameAtt, ncid_out, VarOutId) )
      end do

      call check( nf90_enddef(ncid_out) )

      call get_ndims(var_name,dims)

      if(numDims.eq.1) then
        allocate( lon1d( dims(1) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lon1d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lon1d ) )
      else if(numDims.eq.2) then
        allocate( lon2d( dims(1),dims(2) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lon2d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lon2d ) )
      else if(numDims.eq.3) then
        allocate( lon2d( dims(1),dims(2) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  lon2d ) )
        call check( nf90_put_var(ncid_out, VarOutId, lon2d ) )
      else
        call finalize
        stop "copy_coordinates: UNKNOWN LON FORMAT IN SRC FILE. STOP."
      end if

    
     end if

  end do

  ! Levels
  if(.not.iswrf) then
   do ii = 1, size(lev_names)

    var_name = lev_names(ii)
    status = nf90_inq_varid(ncid_in, var_name, VarInId)

    if(status .eq. nf90_noerr) then 
      
      call check( nf90_redef(ncid_out)  )

      call check( nf90_inquire_variable(ncid_in, VarInId, ndims = numDims,   &
                    dimids = dim_id_var, nAtts = numAtts, xtype = xtype ) )
      ! in WRF coordinates are written in 3D (moving nests accounted)
      ! We do not use moving nests, so converting 3D->2D
      if(numDims.ne.3) then
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:numDims) /), varid = VarOutId ) )
      else
        call check( nf90_def_var(ncid_out, name = var_name, xtype = xtype,     &
                    dimids = (/ dim_id_var(:2) /), varid = VarOutId ) )
      end if

      do ia = 1, numAtts
        call check ( nf90_inq_attname(ncid_in, VarInId, ia, nameAtt) )
        call check ( nf90_copy_att(ncid_in, VarInId, nameAtt, ncid_out, VarOutId) )
      end do

      call check( nf90_enddef(ncid_out) )

      call get_ndims(var_name,dims)

      if(numDims.eq.1) then
        allocate( levels( dims(1) ) )
        call check( nf90_get_var(ncid_in,  VarInId,  levels ) )
        call check( nf90_put_var(ncid_out, VarOutId, levels ) )
      else
        call finalize
        stop "copy_coordinates: UNKNOWN LEV FORMAT IN SRC FILE. STOP."
      end if

    
     end if

  end do
 end if

end subroutine copy_coordinates

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Checks if scaling was applied to the data 
! returns 0 - scale factor is NOT applied
function get_scale(ncid,var_id)
implicit none
    integer,intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=70) :: scale_names(7) = (/"SCALE", "Scale", "_scale", "scale_factor", "Scale_factor", "Slope" , "slope"/)
    logical :: check_scale, scale, offset
    integer ii, xdim,status
    real(kind=real64) :: get_scale

    get_scale = 0

    xdim  = UBOUND(scale_names,1)
    ! if found scale -> stop loop and return the value
    do ii = 1, xdim
        status = nf90_get_att(ncid, var_id, scale_names(ii), get_scale )
        if( status .eq. nf90_noerr ) exit
    end do

end function get_scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Checks if offset was applied to the data 
! returns 0 - offset is NOT applied
function get_offset(ncid,var_id)
implicit none
    integer,intent(in) :: ncid
    integer, intent(in) :: var_id
    character(len=70) :: offset_name(6) = (/"add_offset", "OFFSET", "Offset", "_offset", "Intercept", "intercept"/)
    logical :: check_scale, scale, offset
    integer ii, xdim,status
    real(kind=real64) :: get_offset

    get_offset = 0

    xdim  = UBOUND(offset_name,1)
    ! if found offset -> stop loop and return the value
    do ii = 1, xdim
        status = nf90_get_att(ncid, var_id, offset_name(ii), get_offset )
        if( status .eq. nf90_noerr ) exit
    end do

end function get_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine save_kernel(kernel,sigma,truncate)
    implicit none
    real(kind=real32), intent(in) :: sigma
    real(kind=real32), intent(in), dimension(:,:), allocatable :: kernel
    real(kind=real32), intent(in) :: truncate
    integer             :: ncid_tmp, VarId_tmp, xDimID_tmp, yDimID_tmp
    character(len=80)   :: filename_tmp
    integer             :: unit_num
    integer             :: Mx, My
    real(kind=real32)   :: size

    Mx = ubound(kernel, 1) * 2 + 1
    My = ubound(kernel, 2) * 2 + 1

    size = sigma*truncate*2+1
    
    filename_tmp="_kernel_csi_"//file_in
    
    unit_num = 200 
    open(unit_num, file=trim(filename_tmp), status='unknown')
    close(unit_num, status='delete')
    
    call check (  nf90_create(path=filename_tmp,cmode=or(nf90_clobber,nf90_64bit_offset),ncid=ncid_tmp) )
        call check ( nf90_def_dim(ncid_tmp, "x", Mx, xDimID_tmp) )
        call check ( nf90_def_dim(ncid_tmp, "y", My, yDimID_tmp) ) 
        call check ( nf90_def_var(ncid_tmp, name = "kernel", xtype = NF90_FLOAT, dimids = (/xDimID_tmp,yDimID_tmp/), varid = VarId_tmp) )
        call check ( nf90_put_att(ncid_tmp, NF90_GLOBAL , "Gaussian_Size", size ) )
        call check ( nf90_put_att(ncid_tmp, NF90_GLOBAL , "Sigma", sigma) )
        call check ( nf90_put_att(ncid_tmp, NF90_GLOBAL , "Truncate", truncate) )
        call check ( nf90_put_att(ncid_tmp, NF90_GLOBAL , "Description", "Gaussian_Size=sigma*truncate*2+1") )
    call check (  nf90_enddef(ncid_tmp) )
    call check (  nf90_put_var(ncid_tmp, VarId_tmp, kernel ) )
    call check (  nf90_close(ncid_tmp) )
    
    
end subroutine save_kernel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_stat_3d_32(str,array)
    implicit none
    character(len=*), intent(in) :: str
    real(kind=real32), intent(in), dimension(:,:,:) :: array
    integer :: kk, Mz

    Mz = ubound(array, 3) 

    do kk = 1, Mz
      print*,trim(str),kk,minval(array(:,:,kk),mask=array(:,:,kk).ne.fFillValue),&
      sum(array(:,:,kk),mask=array(:,:,kk).ne.fFillValue)/count(array(:,:,kk).ne.fFillValue),&
      maxval(array(:,:,kk),mask=array(:,:,kk).ne.fFillValue)
    end do
    
end subroutine print_stat_3d_32

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_stat_3d_64(str,array)
    implicit none
    character(len=*), intent(in) :: str
    real(kind=real64), intent(in), dimension(:,:,:) :: array
    integer :: kk, Mz

    Mz = ubound(array, 3) 

    do kk = 1, Mz
      print*,trim(str),kk,minval(array(:,:,kk),mask=array(:,:,kk).ne.dFillValue),&
      sum(array(:,:,kk),mask=array(:,:,kk).ne.dFillValue)/count(array(:,:,kk).ne.dFillValue),&
      maxval(array(:,:,kk),mask=array(:,:,kk).ne.dFillValue)
    end do
    
end subroutine print_stat_3d_64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine finalize
implicit none
  call close_nc(ncid_in)
  call close_nc(ncid_out)
end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine read_namelist
! 
! Reads the namelist. All variables, including name of the namelist
! Are set in module_globals.f90
! 
implicit none
  integer :: unit

  open(newunit=unit,file=trim(namelist_name))
  read(unit,nml=settings)
  read(unit,nml=criteria)
  close(unit)

end subroutine read_namelist

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine check(status)
implicit none
    integer, intent (in) :: status
      
    if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
        stop "Stopped"
    end if
end subroutine check

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine open_src_nc
implicit none

  call check( nf90_open(file_in,  NF90_NOWRITE, ncid_in ) )

end subroutine open_src_nc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine close_nc(ncid)
implicit none
    integer,intent(in) :: ncid
    call check( nf90_close(ncid) )
end subroutine close_nc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module module_io