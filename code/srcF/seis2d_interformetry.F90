program seis2d_interformetry

!----------------------------------------------------------------------------
! This is the main program to apply interformetry to reconstructed wave
! in RTI.
!
! Author: 
!       Yicheng Zhou    Email: zhouyc07@mail.ustc.edu.cn
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------
!
!{ -- declare module used --
use para_mod
use io_mod
use nompi_mod
use media_mod
use src_mod
use interformetry_mod
use interformetry_solver_mod
!} -- end declare module used ---

implicit none
integer :: ntime,ierr
logical :: is_rest
integer :: n
call get_conf_name(fnm_conf)
call interformetry_fnm_init(fnm_conf)
call interformetry_import(fnm_int_conf)

call solver_interformetry_init
call solver_interformetry_get
    print *,'ntimes==',ntimes
    print *,'ntimee==',ntimee
ntime=0


loop_time: do
print *,ntime
  if ((ntime<ntimee) .and.(ntime>ntimes)) then
    print *,'ntime==',ntime
    call update_interformetry(iffo1,iffo2)
    call update_img
  end if
  call wave_interformetry_export(SVpx,SVpz,SVsx,SVsz,GVpx,GVpz,&
                                 GVsx,GVsz,IVpx,IVpz,IVsx,IVsz&
                                 IGVpx,IGVpz,IGVsx,IGVsz,time)
  call img_export(Imgcon,IImgcon)
  ntime=ntime+1
end do loop_time

end program seis2d_interformetry
