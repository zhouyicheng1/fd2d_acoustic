module interformetry_mod
!
!
! This module is used for interformetry
!
! Wei Zhang     Email: zhangwei.zw@gmail.con
! Yicheng Zhou  Email: zhouyc07@mail.ustc.edu.cn
! Copyright(C) 2006 Wei ZHANG          
!
!*********************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $      
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $      
! $Revision: 510 $       
! $LastChangeBy: yicheng zhou $
!
!*********************************************************
!
use constants_mod
use string_mod
use math_mod
use para_mod
use nompi_mod
use grid_mod
use media_mod
use nfseis_mod
implicit none
private

public ::             &
  interformetry_fnm_init,  &
  interformetry_import
character (len=SEIS_STRLEN),public :: fnm_int_conf, fnm_int,pnm_int
integer, public :: iffo1,iffo2,inx1,inx2,inz1,inz2,ntimes,ntimee

!---------------------------------------------------------
contains
!---------------------------------------------------------
subroutine interformetry_fnm_init(filenm)
character (len=*),intent(in) :: filenm
integer fid
fid=1001
open(fid,file=trim(filenm),status="old")
call string_conf(fid,1,'INT_CONF',2,fnm_int_conf)  ! fnm_int_conf --INT_CONF = SeisInterformetry.conf
call string_conf(fid,1,'INT_ROOT',2,pnm_int)       ! int_src --INT_ROOT = ./input
close(fid)
end subroutine interformetry_fnm_init
subroutine interformetry_import(filenm)
character (len=*),intent(in) :: filenm
!character(len =SEIS_STRLEN)::iffo1,iffo2,nxi1,nxi2,nyi1,nyi2,nzi1,nzi2,ntimes,ntimee
integer fid
fid=1001
!open(fid,file=trim(filenm),status="old")
!  call string_conf(fid,1,'INT_CONF',2,fnm_interformetry_conf)  !fnm_src_conf -- SOURCE_CONF = SeisSource.conf
!  call string_conf(fid,1,'INT_ROOT',2,pnm_interformetry)       !pnm_src -- SOURCE_ROOT = ./input
!close(fid)
open(fid,file=trim(filenm),status="old")
call string_conf(fid,1,'time_interformetry_window_start',2,ntimes)
call string_conf(fid,1,'time_interformetry_window_end',2,ntimee)
call string_conf(fid,1,'space_interformetry_window_size1',2,iffo1)
call string_conf(fid,1,'space_interformetry_window_size2',2,iffo2)
call string_conf(fid,1,'x1',2,inx1)
call string_conf(fid,1,'x2',2,inx2)
call string_conf(fid,1,'z1',2,inz1)
call string_conf(fid,1,'z2',2,inz2)
close (fid)
print *, ntimes
print *, ntimee
end subroutine interformetry_import
end module interformetry_mod
