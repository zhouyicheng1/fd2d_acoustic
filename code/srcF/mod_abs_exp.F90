module abs_exp_mod

!----------------------------------------------------------------------------
! This module is used for absorbing outgoing waves based on
! Cerjan and Kosloff's nonreflecting boundary condition
! (Cerjan C., et al.(1985), Geophysics, 50(4):705-708)
! with improvement of taking into account time step value
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------

use constants_mod
use string_mod
use solver_mod, only : P,Vx,Vz
use para_mod
use grid_mod
use nompi_mod
implicit none

private

public ::               &
  abs_exp_init,         &
  abs_exp_destroy,      &
  abs_exp_hook,         &
  abs_exp_momentum,     &
  abs_exp_rest_export,  &
  abs_exp_rest_import

!-----------------------------------------------------------------------------

integer,parameter :: num_blk=4

type ELEM
  logical :: isabs
  integer :: nx1,nx2,nz1,nz2,nx,nz
  integer :: ni1,ni2,nk1,nk2,ni,nk
end type ELEM

type (ELEM), dimension(num_blk) :: W

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: Ex,Ez

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_exp_init(fnm_conf)
  use nompi_mod, only : absnode
  character (len=*),intent(in) :: fnm_conf
  integer :: fid,n,m,i,k
  real(SP),dimension(SEIS_GEO,2) :: Vs

  fid=1001
  abs_number=0
  Vs=0.0_SP
  open(fid,file=trim(fnm_conf),status="old")
  do n=1,SEIS_GEO
  do m=1,2
  if (absnode(n,m)) then
    call string_conf(fid,1,'abs_number',(n-1)*2+m+1,abs_number(n,m))
    ! reset absnode if layer <= 0
    if (abs_number(n,m)<=0) then
      absnode(n,m)=.false.
      cycle
    end if
    call string_conf(fid,1,'abs_velocity',(n-1)*2+m+1,Vs(n,m))
  end if
  end do
  end do
  close(fid)

  ! if layer >0, deset freenode
  if (abs_number(2,2)>0 .and. freenode) then
    freenode=.false.
  end if

  do n=1,num_blk
    W(n)%isabs=.false.
    W(n)%nx1=nx1; W(n)%nx2=nx1-1; W(n)%nx=0
    W(n)%nz1=nz1; W(n)%nz2=nz1-1; W(n)%nz=0
  end do

  allocate(Ex(nx1:nx2)); Ex=1.0_SP
  allocate(Ez(nz1:nz2)); Ez=1.0_SP

  !x1
  n=1
  W(n)%ni=abs_number(1,1);W(n)%ni1=ni1       ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
  W(n)%nk=nk             ;W(n)%nk1=nk1       ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
  do i=W(n)%ni1,W(n)%ni2
    Ex(i)=cal_e(W(n)%ni -(i-W(n)%ni1),Vs(1,1),x(i+1)-x(i),W(n)%ni)
  end do
  if (W(n)%ni>0 .and. W(n)%nk>0) then
    W(n)%isabs=.true.
  end if
  !x2
  n=n+1
  W(n)%ni=abs_number(1,2);W(n)%ni1=ni2-abs_number(1,2)+1;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
  W(n)%nk=nk             ;W(n)%nk1=nk1                  ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
  do i=W(n)%ni1,W(n)%ni2
    Ex(i)=cal_e(i-W(n)%ni1+1         ,Vs(1,2),x(i)-x(i-1),W(n)%ni)
  end do
  if (W(n)%ni>0 .and. W(n)%nk>0) then
    W(n)%isabs=.true.
  end if

  !z1
  n=n+1
  W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1);W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
  W(n)%nk=abs_number(2,1)        ;W(n)%nk1=nk1                ;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
  do k=W(n)%nk1,W(n)%nk2
    Ez(k)=cal_e(W(n)%nk -(k-W(n)%nk1),Vs(2,1),z(k+1)-z(k),W(n)%nk)
  end do
  if (W(n)%ni>0 .and. W(n)%nk>0) then
    W(n)%isabs=.true.
  end if
  !z2
  n=n+1
  W(n)%ni=ni-sum(abs_number(1,:));W(n)%ni1=ni1+abs_number(1,1)  ;W(n)%ni2=W(n)%ni1+W(n)%ni-1; 
  W(n)%nk=abs_number(2,2)        ;W(n)%nk1=nk2-abs_number(2,2)+1;W(n)%nk2=W(n)%nk1+W(n)%nk-1; 
  do k=W(n)%nk1,W(n)%nk2
    Ez(k)=cal_e(k-W(n)%nk1+1         ,Vs(2,2),z(k)-z(k-1),W(n)%nk)
  end do
  if (W(n)%ni>0 .and. W(n)%nk>0) then
    W(n)%isabs=.true.
  end if

  do n=1,num_blk
    if (W(n)%isabs) then
      W(n)%nx1=W(n)%ni1;W(n)%nz1=W(n)%nk1
      W(n)%nx2=W(n)%ni2;W(n)%nz2=W(n)%nk2
    end if
  end do
#ifdef DEBUG
  do i=ni1,ni2
    write(50+myid,*) Ex(i)
  end do
  do k=nk1,nk2
    write(50+myid,*) Ez(k)
  end do
#endif
end subroutine abs_exp_init
subroutine abs_exp_destroy
  return
end subroutine abs_exp_destroy

subroutine abs_exp_momentum
  integer :: n,i,k
  real(SP) :: d

  do n=1,num_blk
  do k=W(n)%nk1,W(n)%nk2
  do i=W(n)%ni1,W(n)%ni2
    d=min(Ex(i),Ez(k))
    Vx (i,k)=Vx (i,k)*d
    Vz (i,k)=Vz (i,k)*d
  end do
  end do
  end do
end subroutine abs_exp_momentum

subroutine abs_exp_hook
  integer :: n,i,k
  real(SP) :: d

  do n=1,num_blk
  do k=W(n)%nk1,W(n)%nk2
  do i=W(n)%ni1,W(n)%ni2
    d=min(Ex(i),Ez(k))
    P(i,k)=P(i,k)*d
  end do
  end do
  end do
end subroutine abs_exp_hook

function cal_e(i,Vs,ah,nb) result(d)
  integer,intent(in) :: i,nb
  real(SP),intent(in) :: Vs,ah
  real(SP) :: d
  real(SP) :: ie
  integer :: m,n
  ie=i
  !Vs=5000.0_SP
  m=(nb*ah)/(Vs*stept)
  d=0.0_SP
  do n=1,m
    d=d+(n*stept*Vs)**2/(nb*ah)**2
  end do
  d=0.99_SP/d*2.0_SP
  d=exp(-d*(ie/nb)**2)
end function cal_e

subroutine abs_exp_rest_export(pnm_rest)
  character (len=*),intent(in) :: pnm_rest
  return
end subroutine abs_exp_rest_export

subroutine abs_exp_rest_import(pnm_rest)
  character (len=*) :: pnm_rest
  return
end subroutine abs_exp_rest_import

end module abs_exp_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
