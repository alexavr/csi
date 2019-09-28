! Description here ...
! 
program csi
  use module_globals
  use module_csi, only: run_computations
  use module_io,  only: print_header, print_footer, read_namelist, open_src_nc, &
                        create_dst_nc, finalize, time_init, time_start, time_end

  implicit none ! religion first

  call get_command_argument(1,file_in)
  if(trim(file_in).EQ."") then
    stop "Usage: ./csi.exe input_file.nc"
  end if
  file_out = "_csi_"//trim(file_in)

  call time_init                      ! Initiate the time count
  call time_start(nb_ticks_initial)   ! Start marker

  call read_namelist                  ! 

  if (debug) call print_header        ! Frint some debug junk 
  
  call open_src_nc                    ! Open the src file

  call create_dst_nc                  ! Create the output file. Also gets IDs 
                                      ! and sizes of dimentions. 
  ! READS VARIABLES, [COMPUTES RHO AND F], COMPUTES GRADIENTS, VG AND VA
  call run_computations ! potentially parallel

  ! DEALLOCATES PERMANENT ARRAYS, CLOSES NC FILES
  call finalize

  call print_footer                   ! Print src filename and elapsed time 

end program csi