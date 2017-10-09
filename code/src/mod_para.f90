

module para_mod

!----------------------------------------------------------------------------
! This module declares parameters variables
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------






































































use constants_mod
use string_mod, only : string_conf
implicit none
!private

public :: para_init,get_conf_name,  &
  para_basic_init,                  &
  inn_i,inn_k,                      &
  out_i,out_k,                      &
  loct_i,loct_k,                    &
  set_cur_time

!-----------------------------------------------------------------------------
integer,public :: NTPI,NTPK,NTPX,NTPZ
integer,public ::                   &
  ni1,ni2,nk1,nk2,ni,nk,            &! without additional stencil points
  nx1,nx2,nz1,nz2,nx,nz,            &! include boundary stencil points
  ngi1,ngi2,ngk1,ngk2,              &! global index without ghost points
  ngx1,ngx2,ngz1,ngz2,              & ! same as above with ghost
  npi1,npi2,npk1,npk2
integer,dimension(SEIS_GEO*2),public :: point_in_this
integer,public :: nt
real(SP),public :: stept
real(SP),public :: cur_time
integer,public :: cur_nt,abs_type
character (len=SEIS_STRLEN),public :: fnm_conf
logical, public ::        &
  MediaPrecise,           &
  SrcSmooth,              &
  SrcForceAverage,        &
  CondFreeAFDA,           &
  CondFreeSIMG,           &
  CondFreeVZERO,          &
  COndFreeV2TH,           &
  CondFreeVEXT,           &
  FreeOnTii,              &
  ResultOnColocatedGrid
!     DataTypeDouble,       &
!     SrcForceAverage,      &
!     MediaPrecise,         &
!     1,          &
!     CondFreeSIMG,         &
!     CondFreeSAFDA,        &
!     CondFreeVAFDA,        &
!     CondFreeVEXT,         &
!     CondFreeVZERO,        &
!     COndFreeV2TH,         &
!     FreeOnTii,            &
!     ResultOnColocatedGrid

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------
subroutine get_conf_name(fnm_input)
  character (*),intent(out) :: fnm_input

  fnm_input="SeisFD2D.conf"
end subroutine get_conf_name

!---------------------------------------------------------------------------
subroutine para_basic_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid

  fid=1002
  open(fid,file=trim(fnm_conf),status="old")
  !About media
  call string_conf(fid,1,'MediaPrecise' ,2,MediaPrecise)
  !About source
  call string_conf(fid,1,'SrcSmooth' ,2,SrcSmooth)
  call string_conf(fid,1,'SrcForceAverage' ,2,SrcForceAverage)
  !About free condition (AFDA)
  call string_conf(fid,1,'CondFreeAFDA',2,CondFreeAFDA)
  !About free condition on stress
  call string_conf(fid,1,'CondFreeSIMG' ,2,CondFreeSIMG)
  !About free condition on velocity
  call string_conf(fid,1,'CondFreeVZERO',2,CondFreeVZERO)
  call string_conf(fid,1,'CondFreeV2TH',2,CondFreeV2TH)
  call string_conf(fid,1,'CondFreeVEXT',2,CondFreeVEXT)
  !About free-surface position
  call string_conf(fid,1,'FreeOnTii',2,FreeOnTii)
  !About output result
  call string_conf(fid,1,'ResultOnColocatedGrid',2,ResultOnColocatedGrid)
  print *, "MediaPrecise            is ", MediaPrecise
  print *, "SrcSmooth               is ", SrcSmooth
  print *, "SrcForceAverage         is ", SrcForceAverage
  print *, "CondFreeAFDA            is ", CondFreeAFDA
  print *, "CondFreeSIMG            is ", CondFreeSIMG
  print *, "CondFreeVZERO           is ", CondFreeVZERO
  print *, "CondFreeV2TH            is ", CondFreeV2TH
  print *, "CondFreeVEXT            is ", CondFreeVEXT
  print *, "FreeOnTii               is ", FreeOnTii
  print *, "ResultOnColocatedGrid   is ", ResultOnColocatedGrid
  close(fid)
end subroutine

subroutine para_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid
  abs_type=0

  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
    call string_conf(fid,1,'ni',2,ni)
    call string_conf(fid,1,'nk',2,nk)
    call string_conf(fid,1,'nt',2,nt)
    nx1=1; ni1=nx1+2; ni2=ni1+ni-1; nx2=ni2+2
    nz1=1; nk1=nz1+2; nk2=nk1+nk-1; nz2=nk2+2
    nx=nx2-nx1+1; nz=nz2-nz1+1

    call string_conf(fid,1,'stept',2,stept)
    call string_conf(fid,1,'abs_type',2,abs_type);

  close(fid)
  cur_time=0.0
  cur_nt  = 0
end subroutine para_init

subroutine reset_nt(ntime)
  integer,intent(in) :: ntime
  nt = ntime
  print *, 'reset nt to ',nt,' step'
end subroutine reset_nt

subroutine set_cur_time(ntime,incr)
  integer,intent(in) :: ntime
  real(SP),intent(in) :: incr
  cur_nt = ntime
  cur_time = (ntime+incr)*stept
end subroutine set_cur_time

function inn_i(oi) result(ii)
  integer,intent(in) :: oi
  integer :: ii
  ii=oi+ni1-1
end function inn_i
function inn_k(ok) result(kk)
  integer,intent(in) :: ok
  integer :: kk
  kk=ok+nk1-1
end function inn_k

function out_i(ii) result(oi)
  integer,intent(in) :: ii
  integer :: oi
  oi=ii-ni1+1
end function out_i
function out_k(kk) result(ok)
  integer,intent(in) :: kk
  integer :: ok
  ok=kk-nk1+1
end function out_k
function loct_i(n) result(i)
  integer,intent(in) :: n
  integer :: i
  i=n+nx1-1
end function loct_i
function loct_k(n) result(k)
  integer,intent(in) :: n
  integer :: k
  k=n+nz1-1
end function loct_k

end module para_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
