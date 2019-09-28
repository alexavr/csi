module module_csi
use module_globals
use module_io, only:  get_nc_var, wrf_user_unstagger, put_var_time, save_kernel, &
                      print_stat, finalize, time_start, time_end

implicit none

private

public :: run_computations 

contains 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine run_computations ! potentially parallel
  integer :: itime, nb_ticks_run_computations 
  real(kind=real32) :: elapsed_time_run_computations
  logical :: first_time
  real(kind=real64),dimension(:,:,:),allocatable :: dudx, dudy, dudz,   &
                                                    dvdx, dvdy, dvdz,   &
                                                    dwdx, dwdy, dwdz

  ! calculate gradients (center difference)
  allocate( dudx ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dudy ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dudz ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dvdx ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dvdy ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dvdz ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dwdx ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dwdy ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( dwdz ( xDimSIZE, yDimSIZE, zDimSIZE ) )

  first_time = .TRUE.

  do itime = 1, tDimSIZE

    if (debug) call time_start(nb_ticks_run_computations)
    
    ! this subs read data and make stress tensor derivatives (dudx, etc...) 
    ! based on the grid. taking care of staggering grid in WRF.
    if(iswrf) then

      call prep_data_wrf(itime,first_time,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) 
      
    else
      
      call prep_data_rean(itime,first_time,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)
    
    end if

    if (q_criterion .OR. delta_criterion .OR. l2_criterion) then
      call compute_eurlean_criteria(itime,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) 
    end if

    if (debug) then
      call time_end(nb_ticks_run_computations,elapsed_time_run_computations)
      print'(a,i3,a,i5,a,f7.2,a)',"run_computations: time step ",itime," of ", tDimSIZE, " took ", elapsed_time_run_computations, " sec"
    end if

    first_time = .FALSE.

  end do

  deallocate(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

end subroutine run_computations

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

subroutine prep_data_wrf(itime,first_time,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) ! potentially parallel
use gaussian_filter, only: gaussian_kernel, convolve
  integer, intent(in) :: itime 
  logical, intent(in) :: first_time 
  real(kind=real64),dimension(:,:,:),intent(inout) :: dudx,dudy,dudz,&
                                                      dvdx,dvdy,dvdz,&
                                                      dwdx,dwdy,dwdz 
  real(kind=real32),dimension(:,:)  ,allocatable :: MAPX,MAPY
  real(kind=real64),dimension(:,:)  ,allocatable :: RDX,RDY 
  real(kind=real64),dimension(:,:,:),allocatable :: RDZ 
  real(kind=real32),dimension(:,:,:),allocatable :: U,V,W,PH 
  real(kind=real32),dimension(:,:,:),allocatable :: tmp3d1,tmp3d2,tmp3d3
  real(kind=real32),allocatable :: kernel(:,:)
  real(kind=real32) :: RDX1D,RDY1D
  integer :: ii,jj,kk
  integer :: dimIDs_m(4)

  dimIDs_m = (/ xDimIDout, yDimIDout, zDimIDout, tDimIDout /)

  dudx = dFillValue
  dudy = dFillValue
  dudz = dFillValue
  dvdx = dFillValue
  dvdy = dFillValue
  dvdz = dFillValue
  dwdx = dFillValue
  dwdy = dFillValue
  dwdz = dFillValue

  allocate( RDX  ( xDimSIZE, yDimSIZE ) )
  allocate( RDY  ( xDimSIZE, yDimSIZE ) )

  if(first_time) then
    allocate( MAPX ( xDimSIZE, yDimSIZE ) )
    allocate( MAPY ( xDimSIZE, yDimSIZE ) )

    call get_nc_var(ncid_in,"RDX",itime,RDX1D)
    call get_nc_var(ncid_in,"RDY",itime,RDY1D)
    call get_nc_var(ncid_in,"MAPFAC_MX",itime,MAPX)
    call get_nc_var(ncid_in,"MAPFAC_MY",itime,MAPY)

    RDX = RDX1D*1.d0/MAPX ! 1.d0/DBLE(DX*MAPX)
    RDY = RDY1D*1.d0/MAPY ! 1.d0/DBLE(DY*MAPY)

    deallocate( MAPX,MAPY )
  end if

  ! getting U (plus unstaggering)
  allocate( U ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( tmp3d1 ( xDimSIZE_stag, yDimSIZE, zDimSIZE ) )
  call get_nc_var(ncid_in, "U", itime, tmp3d1)
  call wrf_user_unstagger(tmp3d1,U)      
  deallocate( tmp3d1 )

  ! getting V (plus unstaggering)
  allocate( V ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( tmp3d1 ( xDimSIZE, yDimSIZE_stag, zDimSIZE ) )
  call get_nc_var(ncid_in, "V", itime, tmp3d1)
  call wrf_user_unstagger(tmp3d1,V)      
  deallocate( tmp3d1 )

  ! getting W (plus unstaggering)
  allocate( W ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( tmp3d1 ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  call get_nc_var(ncid_in, "W", itime, tmp3d1)
  call wrf_user_unstagger(tmp3d1,W)      
  deallocate( tmp3d1 )

  if(debug_output) then
    call put_var_time(ncid_out,U,"U","U","m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,V,"V","V","m s-1",itime,dimIDs_m)
    call put_var_time(ncid_out,W,"W","W","m s-1",itime,dimIDs_m)
  end if

  if(smooth) then
    
    call gaussian_kernel(smooth_sigma, kernel, smooth_truncate) ! GAVR: Not efficient! Kernel shouldn't be calculated every timestep
    
    allocate( tmp3d1 ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    allocate( tmp3d2 ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    allocate( tmp3d3 ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    
    do kk = 1, zDimSIZE
      call convolve(U(:,:,kk), kernel, tmp3d1(:,:,kk))
      call convolve(V(:,:,kk), kernel, tmp3d2(:,:,kk))
      call convolve(W(:,:,kk), kernel, tmp3d3(:,:,kk))
    end do

    U = tmp3d1
    V = tmp3d2
    W = tmp3d3

    if(debug_output) then
      call put_var_time(ncid_out,U,"GU","U","m s-1",itime,dimIDs_m)
      call put_var_time(ncid_out,V,"GV","V","m s-1",itime,dimIDs_m)
      call put_var_time(ncid_out,W,"GW","W","m s-1",itime,dimIDs_m)
      if(first_time) call save_kernel(kernel,smooth_sigma,smooth_truncate)
    end if

    deallocate(kernel)
    deallocate(tmp3d1,tmp3d2,tmp3d3)

  end if

  ! print*,sum(U(:,:,1))/(xDimSIZE*yDimSIZE),&
  !        sum(V(:,:,1))/(xDimSIZE*yDimSIZE),&
  !        sum(W(:,:,1))/(xDimSIZE*yDimSIZE)

  ! inverse vertical level steps (1/m)
  allocate( RDZ ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  allocate( tmp3d1 ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  allocate( tmp3d2 ( xDimSIZE, yDimSIZE, zDimSIZE_stag ) )
  call get_nc_var(ncid_in, "PH",  itime, tmp3d1 )
  call get_nc_var(ncid_in, "PHB", itime, tmp3d2 )
  tmp3d2 = (tmp3d1 + tmp3d2)/g
  do kk=1, zDimSIZE
      RDZ(:,:,kk) = 1.d0/DBLE(tmp3d2(:,:,kk+1) - tmp3d2(:,:,kk))
  end do

  if(output_height) then
    allocate( PH  ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    call wrf_user_unstagger(tmp3d2,PH)
    call put_var_time(ncid_out,PH,"Z","PH/g","m",itime,dimIDs_m)
    deallocate( PH )
  end if

  deallocate( tmp3d1, tmp3d2 )

  do ii=1, xDimSIZE-1
    do jj=1, yDimSIZE-1
      do kk=1, zDimSIZE-1

        dudx(ii,jj,kk) = RDX(ii,jj)    * DBLE(U(ii+1,jj,kk) - U(ii,jj,kk))
        dudy(ii,jj,kk) = RDY(ii,jj)    * DBLE(U(ii,jj+1,kk) - U(ii,jj,kk))
        dudz(ii,jj,kk) = RDZ(ii,jj,kk) * DBLE(U(ii,jj,kk+1) - U(ii,jj,kk)) 

        dvdx(ii,jj,kk) = RDX(ii,jj)    * DBLE(V(ii+1,jj,kk) - V(ii,jj,kk))
        dvdy(ii,jj,kk) = RDY(ii,jj)    * DBLE(V(ii,jj+1,kk) - V(ii,jj,kk))
        dvdz(ii,jj,kk) = RDZ(ii,jj,kk) * DBLE(V(ii,jj,kk+1) - V(ii,jj,kk)) 

        dwdx(ii,jj,kk) = RDX(ii,jj)    * DBLE(W(ii+1,jj,kk) - W(ii,jj,kk))
        dwdy(ii,jj,kk) = RDY(ii,jj)    * DBLE(W(ii,jj+1,kk) - W(ii,jj,kk))
        dwdz(ii,jj,kk) = RDZ(ii,jj,kk) * DBLE(W(ii,jj,kk+1) - W(ii,jj,kk)) 

      end do
    end do
  end do

  deallocate(U,V,W)
  deallocate(RDX,RDY,RDZ)

  if(debug_output) then
    call put_var_time(ncid_out,dudx,"dudx","dudx","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dudy,"dudy","dudy","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dudz,"dudz","dudz","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dvdx,"dvdx","dvdx","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dvdy,"dvdy","dvdy","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dvdz,"dvdz","dvdz","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dwdx,"dwdx","dwdx","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dwdy,"dwdy","dwdy","m s-2",itime,dimIDs_m)
    call put_var_time(ncid_out,dwdz,"dwdz","dwdz","m s-2",itime,dimIDs_m)
  end if


end subroutine prep_data_wrf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! Assume reanalyze is on the p-levels 
subroutine prep_data_rean(itime,first_time,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) 
  integer, intent(in) :: itime 
  logical, intent(in) :: first_time 
  real(kind=real64),dimension(:,:,:),intent(inout) :: dudx,dudy,dudz,&
                                                      dvdx,dvdy,dvdz,&
                                                      dwdx,dwdy,dwdz 

end subroutine prep_data_rean

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
subroutine compute_eurlean_criteria(itime, dudx, dudy, dudz, dvdx, &
                                    dvdy, dvdz, dwdx, dwdy, dwdz) 
  integer, intent(in) :: itime 
  real(kind=real64),dimension(:,:,:),intent(in)  :: dudx,dudy,dudz,  &
                                                    dvdx,dvdy,dvdz,  &
                                                    dwdx,dwdy,dwdz 
  real(kind=real32),dimension(:,:,:),allocatable :: Q, D, L2 
  real(kind=real64),dimension(:,:,:),allocatable :: dQ, dD, dL2
  real(kind=real64),dimension(3,3) :: S, SS, O, OO, A, Y, C
  real(kind=real64),dimension(3)   :: B
  real(kind=real64) :: kd, kdx1, kdx2, kdx3, ka, kb, kc, kq, kp, kt, x1, x2, x3
  real(kind=real64) :: R, norm_O, norm_S
  real(kind=real64) :: normalization
  integer :: dimIDs_m(4)
  integer :: ii,jj,kk

  dimIDs_m = (/ xDimIDout, yDimIDout, zDimIDout, tDimIDout /)

  if(q_criterion    ) allocate( dQ ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  if(delta_criterion) allocate( dD ( xDimSIZE, yDimSIZE, zDimSIZE ) )
  if(l2_criterion   ) allocate( dL2( xDimSIZE, yDimSIZE, zDimSIZE ) )

  if(q_criterion    ) dQ  = dFillValue ! [s-2]
  if(delta_criterion) dD  = dFillValue ! [s-6]
  if(l2_criterion   ) dL2 = dFillValue ! [s-4]

  do kk = 1, zDimSIZE-1
  do ii = 1, xDimSIZE-1
  do jj = 1, yDimSIZE-1

      ! R = -det(D) [velocity gradient tensor]
      R = - ( dudx(ii,jj,kk)*dvdy(ii,jj,kk)*dwdz(ii,jj,kk) + &
              dudy(ii,jj,kk)*dvdz(ii,jj,kk)*dwdx(ii,jj,kk) + &
              dudz(ii,jj,kk)*dvdx(ii,jj,kk)*dwdy(ii,jj,kk) - &
              dudz(ii,jj,kk)*dvdy(ii,jj,kk)*dwdx(ii,jj,kk) - &
              dudy(ii,jj,kk)*dvdx(ii,jj,kk)*dwdz(ii,jj,kk) - &
              dudx(ii,jj,kk)*dvdz(ii,jj,kk)*dwdy(ii,jj,kk) )

      ! S - rate-of-strain tensor
      S(1,1) = dudx(ii,jj,kk); S(1,2) = 0.5d0*(dudy(ii,jj,kk)+dvdx(ii,jj,kk)); S(1,3) = 0.5d0*(dudz(ii,jj,kk)+dwdx(ii,jj,kk)); 
      S(2,1) = S(1,2);         S(2,2) = dvdy(ii,jj,kk);                        S(2,3) = 0.5d0*(dvdz(ii,jj,kk)+dwdy(ii,jj,kk)); 
      S(3,1) = S(1,3);         S(3,2) = S(2,3);                                S(3,3) = dwdz(ii,jj,kk); 
      
      ! O - vorticity tensor
      O(1,1) = 0d0    ; O(1,2) = 0.5d0*(dudy(ii,jj,kk)-dvdx(ii,jj,kk)); O(1,3) = 0.5d0*(dudz(ii,jj,kk)-dwdx(ii,jj,kk)); 
      O(2,1) = -O(1,2); O(2,2) = 0d0;                                   O(2,3) = 0.5d0*(dvdz(ii,jj,kk)-dwdy(ii,jj,kk)); 
      O(3,1) = -O(1,3); O(3,2) = -O(2,3);                               O(3,3) = 0d0; 

      ! norms O and S
      norm_O = sum(O**2)
      norm_S = sum(S**2)

      ! Q and D criteria
      if(q_criterion    ) dQ(ii,jj,kk) = 0.5d0*( norm_O - norm_S ) 
      if(delta_criterion) dD(ii,jj,kk) = dQ(ii,jj,kk)**3/27.d0 + R**2/4.d0 

      if(l2_criterion) then

        ! SS tensor
        SS(1,1) = S(1,1)**2     + S(1,2)*S(2,1) + S(1,3)*S(3,1);  
        SS(1,2) = S(1,1)*S(1,2) + S(1,2)*S(2,2) + S(1,3)*S(3,2);  
        SS(1,3) = S(1,1)*S(1,3) + S(1,2)*S(2,3) + S(1,3)*S(3,3); 
        SS(2,1) = S(2,1)*S(1,1) + S(2,2)*S(2,1) + S(2,3)*S(3,1);  
        SS(2,2) = S(2,1)*S(1,2) + S(2,2)**2     + S(2,3)*S(3,2);  
        SS(2,3) = S(2,1)*S(1,3) + S(2,2)*S(2,3) + S(2,3)*S(3,3);  
        SS(3,1) = S(3,1)*S(1,1) + S(3,2)*S(2,1) + S(3,3)*S(3,1);  
        SS(3,2) = S(3,1)*S(1,2) + S(3,2)*S(2,2) + S(3,3)*S(3,2);  
        SS(3,3) = S(3,1)*S(1,3) + S(3,2)*S(2,3) + S(3,3)**2;  

        ! OO tensor
        OO(1,1) = O(1,1)**2     + O(1,2)*O(2,1) + O(1,3)*O(3,1);  
        OO(1,2) = O(1,1)*O(1,2) + O(1,2)*O(2,2) + O(1,3)*O(3,2);  
        OO(1,3) = O(1,1)*O(1,3) + O(1,2)*O(2,3) + O(1,3)*O(3,3); 
        OO(2,1) = O(2,1)*O(1,1) + O(2,2)*O(2,1) + O(2,3)*O(3,1);  
        OO(2,2) = O(2,1)*O(1,2) + O(2,2)**2     + O(2,3)*O(3,2);  
        OO(2,3) = O(2,1)*O(1,3) + O(2,2)*O(2,3) + O(2,3)*O(3,3);  
        OO(3,1) = O(3,1)*O(1,1) + O(3,2)*O(2,1) + O(3,3)*O(3,1);  
        OO(3,2) = O(3,1)*O(1,2) + O(3,2)*O(2,2) + O(3,3)*O(3,2);  
        OO(3,3) = O(3,1)*O(1,3) + O(3,2)*O(2,3) + O(3,3)**2;  

        A = OO + SS

        ! Kirilov method [A. Dernov]
        Y(1,1) = A(1,1)
        Y(1,2) = A(2,1)
        Y(1,3) = A(3,1)
        Y(2,1) = Y(1,1)*A(1,1) + Y(1,2)*A(1,2) + Y(1,3)*A(1,3)
        Y(2,2) = Y(1,1)*A(2,1) + Y(1,2)*A(2,2) + Y(1,3)*A(2,3)
        Y(2,3) = Y(1,1)*A(3,1) + Y(1,2)*A(3,2) + Y(1,3)*A(3,3)
        Y(3,1) = Y(2,1)*A(1,1) + Y(2,2)*A(1,2) + Y(2,3)*A(1,3)
        Y(3,2) = Y(2,1)*A(2,1) + Y(2,2)*A(2,2) + Y(2,3)*A(2,3)
        Y(3,3) = Y(2,1)*A(3,1) + Y(2,2)*A(3,2) + Y(2,3)*A(3,3)
        ! matrix C - SLAU
        C(1,1) = Y(2,1); C(1,2) = Y(1,1); C(1,3) = 1.d0; 
        C(2,1) = Y(2,2); C(2,2) = Y(1,2); C(2,3) = 0.d0; 
        C(3,1) = Y(2,3); C(3,2) = Y(1,3); C(3,3) = 0.d0; 
        ! B
        B(1) = -Y(3,1); B(2) = -Y(3,2); B(3) = -Y(3,3);
        ! Kramer method 
        kd   = C(1,1)*C(2,2)*C(3,3) + C(1,2)*C(2,3)*C(3,1) + C(1,3)*C(2,1)*C(3,2) - C(1,3)*C(2,2)*C(3,1) - C(1,1)*C(2,3)*C(3,2) - C(1,2)*C(2,1)*C(3,3) 
        kdx1 =   B(1)*C(2,2)*C(3,3) + C(1,2)*C(2,3)*B(3)   + C(1,3)*B(2)*C(3,2)   - C(1,3)*C(2,2)*B(3)   -   B(1)*C(2,3)*C(3,2) - C(1,2)*B(2)*C(3,3) 
        kdx2 = C(1,1)*B(2)*C(3,3)   + B(1)  *C(2,3)*C(3,1) + C(1,3)*C(2,1)*B(3)   - C(1,3)*B(2)*C(3,1)   - C(1,1)*C(2,3)*B(3)   - B(1)*C(2,1)*C(3,3) 
        kdx3 = C(1,1)*C(2,2)*B(3)   + C(1,2)*B(2)  *C(3,1) + B(1)*C(2,1)*C(3,2)   - B(1)*C(2,2)*C(3,1)   - C(1,1)*B(2)*C(3,2)   - C(1,2)*C(2,1)*B(3) 

        ka = kdx1/kd
        kb = kdx2/kd
        kc = kdx3/kd

        kq = ( ka**2 - 3.d0*kb ) / 9.D0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kp = ( 2.d0*ka**3 - 9.d0*ka*kb + 27.d0*kc) / 54.d0
        kt = acos( kp/dsqrt(kq**3) ) / 3.d0

        x1 = -2.d0*dsqrt(kq)*cos(kt)                  - ka/3.d0
        x2 = -2.d0*dsqrt(kq)*cos(kt + (2.d0*pi/3.d0)) - ka/3.d0
        x3 = -2.d0*dsqrt(kq)*cos(kt - (2.d0*pi/3.d0)) - ka/3.d0

        dL2(ii,jj,kk) = maxval( (/x1,x2,x3/),mask=(/x1,x2,x3/).NE.maxval((/x1,x2,x3/)) ) ! [s-4]

      end if

  end do
  end do
  end do

  ! output
  if(q_criterion) then

    allocate( Q ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    Q = fFillValue
    
    ! call print_stat("dQ: ",dQ)

    ! Q_norma = sum(dQ,mask=dQ.ne.dFillValue)/count(dQ.ne.dFillValue)
    normalization = maxval(dQ,mask=dQ.ne.dFillValue)
    
    where(dQ.ne.dFillValue) Q = dQ/normalization 

    call put_var_time(ncid_out,Q,"Q","Q criterion","",itime,dimIDs_m)

    ! call print_stat("Q: ",Q)

    deallocate( dQ, Q )

  end if

  if(delta_criterion) then

    allocate( D ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    D = fFillValue

    normalization = maxval(dD,mask=dD.ne.dFillValue)

    where(dD.ne.dFillValue) D = dD/normalization 
    call put_var_time(ncid_out,D,"D","Delta criterion","",itime,dimIDs_m)

    deallocate( dD, D )

  end if

  if(l2_criterion) then

    allocate( L2 ( xDimSIZE, yDimSIZE, zDimSIZE ) )
    L2 = fFillValue

    normalization = maxval(dL2,mask=dL2.ne.dFillValue)

    where(dL2.ne.dFillValue) L2 = dL2/normalization 
    call put_var_time(ncid_out,L2,"L2","L2 criterion","",itime,dimIDs_m)

    deallocate( dL2, L2 )
    
  end if

end subroutine compute_eurlean_criteria

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

end module module_csi