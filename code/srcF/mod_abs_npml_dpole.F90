module abs_dpcfs_mod

!----------------------------------------------------------------------------
! This module is used for absorbing outgoing waves based on
! unsplit-field CFS ADE-PML (auxiliary differential equation
! complex frequncy shifted PML)
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Haike FENG      Email: fenghk@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------

#if defined FourthORDER
#include "mod_staggerfd4.h"
#elif defined NONUNIFORMGRID
#include "mod_staggerfd4_nonuni.h"
#else
#include "mod_staggerfd2.h"
#endif

use constants_mod
use string_mod, only : string_conf
use para_mod
use nfseis_mod
use solver_mod, only : &
    P,Vx,Vz
use media_mod
use grid_mod
use nompi_mod
implicit none

private
public ::                 &
  abs_dpcfs_init,         &
  abs_dpcfs_destroy,      &
  abs_dpcfs_hook,         &
  abs_dpcfs_momentum,     &
  abs_dpcfs_rest_export,  &
  abs_dpcfs_rest_import

!-----------------------------------------------------------------------------
DEFFDWET 

!integer,parameter :: SP=kind(1.0)
!real(SP),parameter :: CONSPD=2.0_SP
!real(SP),parameter :: CONSPB=2.0_SP
!real(SP),parameter :: CONSPA=1.0_SP

real(SP),dimension(:,:),allocatable ::   &
  Vx11a,Vz11a,P11a, &
  Vx11b,Vz11b,P11b
real(SP),dimension(:,:),allocatable ::   &
  Vx31a,Vz31a,P31a, &
  Vx31b,Vz31b,P31b

real(SP),dimension(:,:),allocatable ::   &
  Vx12a,Vz12a,P12a, &
  Vx12b,Vz12b,P12b
real(SP),dimension(:,:),allocatable ::   &
  Vx32a,Vz32a,P32a, &
  Vx32b,Vz32b,P32b

logical :: isx1,isx2,isz1,isz2

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: A1x_regu,A1z_regu
real(SP),dimension(:),allocatable :: A1x_half,A1z_half
real(SP),dimension(:),allocatable :: B1x_regu,B1z_regu
real(SP),dimension(:),allocatable :: B1x_half,B1z_half
real(SP),dimension(:),allocatable :: D1x_regu,D1z_regu
real(SP),dimension(:),allocatable :: D1x_half,D1z_half
real(SP),dimension(:),allocatable :: A2x_regu,A2z_regu
real(SP),dimension(:),allocatable :: A2x_half,A2z_half
real(SP),dimension(:),allocatable :: B2x_regu,B2z_regu
real(SP),dimension(:),allocatable :: B2x_half,B2z_half
real(SP),dimension(:),allocatable :: D2x_regu,D2z_regu
real(SP),dimension(:),allocatable :: D2x_half,D2z_half

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_dpcfs_init(fnm_conf)
  use nompi_mod, only : absnode,freenode
  character (len=*) :: fnm_conf
  integer fid,n,m,i,k,ierr,npt
  real(SP),dimension(SEIS_GEO,2) :: Vs,fc,Rpp,b1max,d1max,a1max,b2max,d2max,a2max
  real(SP),dimension(SEIS_GEO,2) :: scaled1,scaled2,scalea1,scalea2
  real(SP) :: x0,z0,L0,Lx,Lz
  integer a1type,a2type,d1type,d2type
  real(SP) :: CONSPD1,CONSPD2,CONSPB1,CONSPB2,CONSPA1,CONSPA2

  fid=1001
  abs_number=0
  Vs=0.0_SP
  Rpp=0.0_SP
  d1max=0.0_SP
  a1max=0.0_SP
  b1max=0.0_SP
  d2max=0.0_SP
  a2max=0.0_SP
  b2max=0.0_SP
  scaled1=0.0_SP
  scaled2=0.0_SP
  scalea1=0.0_SP
  CONSPD1=0.0_SP
  CONSPD2=0.0_SP
  CONSPB1=0.0_SP
  CONSPB2=0.0_SP
  CONSPA1 =1.0_SP
  CONSPA2 =1.0_SP

  isx1=.false.
  isx2=.false.
  isz1=.false.
  isz2=.false.

  open(fid,file=trim(fnm_conf),status="old")
  do n=1,SEIS_GEO
  do m=1,2
    if (absnode(n,m)) then
      call string_conf(fid,1,'abs_number',(n-1)*2+m+1,abs_number(n,m))
      !=== reset absnode if layer <= 0 ===
      if (abs_number(n,m)<=0) then
        absnode(n,m)=.false.
        cycle
      end if
      call string_conf(fid,1,'abs_velocity',(n-1)*2+m+1,Vs(n,m))
      call string_conf(fid,1,'CFS_b1max',(n-1)*2+m+1,b1max(n,m))
      call string_conf(fid,1,'CFS_a1max',(n-1)*2+m+1,a1max(n,m))
      call string_conf(fid,1,'CFS_a2max',(n-1)*2+m+1,a2max(n,m))
      call string_conf(fid,1,'CFS_b2max',(n-1)*2+m+1,b2max(n,m))
      call string_conf(fid,1,'CFS_scaled1',(n-1)*2+m+1,scaled1(n,m))
      call string_conf(fid,1,'CFS_scaled2',(n-1)*2+m+1,scaled2(n,m))
      call string_conf(fid,1,'CFS_scalea1',(n-1)*2+m+1,scalea1(n,m))
      call string_conf(fid,1,'CFS_scalea2',(n-1)*2+m+1,scalea2(n,m))
      if (b1max(n,m)<1.0) then
        print *, "n,m,CFS_b1max(n,m)=",n,m,b1max(n,m)
        call swmpi_except("CFS_b1max should be large or equal to 1")
      end if
      Rpp(n,m)=cal_pml_R(abs_number(n,m))
    end if
  end do
  end do
  call string_conf(fid,1,'CONSPD1',2,CONSPD1)
  call string_conf(fid,1,'CONSPD2',2,CONSPD2)
  call string_conf(fid,1,'CONSPB1',2,CONSPB1)
  call string_conf(fid,1,'CONSPB2',2,CONSPB2)
  call string_conf(fid,1,'CONSPA1',2,CONSPA1)
  call string_conf(fid,1,'CONSPA2',2,CONSPA2)
  call string_conf(fid,1,'A1Type',  2,  a1type)
  call string_conf(fid,1,'A2Type',  2,  a2type)
  call string_conf(fid,1,'D1Type',  2,  d1type)
  call string_conf(fid,1,'D2Type',  2,  d2type)
  close(fid)

  !=== if layer >0, unset freenode ==
  if (abs_number(2,2)>0 .and. freenode) then
    freenode=.false.
  end if

  allocate(A1x_regu(nx1:nx2)); A1x_regu=0.0_SP
  allocate(B1x_regu(nx1:nx2)); B1x_regu=1.0_SP
  allocate(D1x_regu(nx1:nx2)); D1x_regu=0.0_SP
  allocate(A1z_regu(nz1:nz2)); A1z_regu=0.0_SP
  allocate(B1z_regu(nz1:nz2)); B1z_regu=1.0_SP
  allocate(D1z_regu(nz1:nz2)); D1z_regu=0.0_SP
  allocate(A1x_half(nx1:nx2)); A1x_half=0.0_SP
  allocate(B1x_half(nx1:nx2)); B1x_half=1.0_SP
  allocate(D1x_half(nx1:nx2)); D1x_half=0.0_SP
  allocate(A1z_half(nz1:nz2)); A1z_half=0.0_SP
  allocate(B1z_half(nz1:nz2)); B1z_half=1.0_SP
  allocate(D1z_half(nz1:nz2)); D1z_half=0.0_SP
  allocate(A2x_regu(nx1:nx2)); A2x_regu=0.0_SP
  allocate(B2x_regu(nx1:nx2)); B2x_regu=1.0_SP
  allocate(D2x_regu(nx1:nx2)); D2x_regu=0.0_SP
  allocate(A2z_regu(nz1:nz2)); A2z_regu=0.0_SP
  allocate(B2z_regu(nz1:nz2)); B2z_regu=1.0_SP
  allocate(D2z_regu(nz1:nz2)); D2z_regu=0.0_SP
  allocate(A2x_half(nx1:nx2)); A2x_half=0.0_SP
  allocate(B2x_half(nx1:nx2)); B2x_half=1.0_SP
  allocate(D2x_half(nx1:nx2)); D2x_half=0.0_SP
  allocate(A2z_half(nz1:nz2)); A2z_half=0.0_SP
  allocate(B2z_half(nz1:nz2)); B2z_half=1.0_SP
  allocate(D2z_half(nz1:nz2)); D2z_half=0.0_SP

  !== the PML starts from regular grid, where normal stress defined ==

  ! -------------------  x layer --------------------------
  ! PML layer toward -x direction, starting from half point where Vx located
  if (abs_number(1,1)>0) then
    isx1=.true.
    x0=(x(abs_number(1,1)+ni1+1)+x(abs_number(1,1)+ni1))/2.0_SP
    L0=x0-x(ni1)
    d1max(1,1)=cal_pml_dmax(L0,Vs(1,1),Rpp(1,1),b1max(1,1),CONSPD1,scaled1(1,1),a2max(1,1),d1type)
    d2max(1,1)=cal_pml_dmax(L0,Vs(1,1),Rpp(1,1),b2max(1,1),CONSPD2,scaled2(1,1),a1max(1,1),d2type)
    !  dmax(1,1)=0.0_SP
    do i=ni1,abs_number(1,1)+ni1
      ! regular grid point
      Lx=x0-x(i)
      D1x_regu(i)=cal_pml_d(Lx,L0,d1max(1,1),CONSPD1)
      D2x_regu(i)=cal_pml_d(Lx,L0,d2max(1,1),CONSPD2)
      B1x_regu(i)=cal_pml_b(Lx,L0,b1max(1,1),CONSPB1)
      B2x_regu(i)=cal_pml_b(Lx,L0,b2max(1,1),CONSPB2)
      A1x_regu(i)=cal_pml_a(a1type,Lx,L0, d2max(1,1),a1max(1,1),scalea1(1,1),CONSPA1)
      !A1x_regu(i)=cal_pml_a(atype,Lx,L0,D2x_regu(i),a1max(1,1),scalea1(1,1),CONSPA1)
      A2x_regu(i)=cal_pml_a(a2type,Lx,L0, d1max(1,1),a2max(1,1),scalea2(1,1),CONSPA2)
      ! half grid point
      Lx=x0-(x(i)+x(i+1))/2.0_SP
      D1x_half(i)=cal_pml_d(Lx,L0,d1max(1,1),CONSPD1)
      D2x_half(i)=cal_pml_d(Lx,L0,d2max(1,1),CONSPD2)
      B1x_half(i)=cal_pml_b(Lx,L0,b1max(1,1),CONSPB1)
      B2x_half(i)=cal_pml_b(Lx,L0,b2max(1,1),CONSPB2)
      A1x_half(i)=cal_pml_a(a1type,Lx,L0, d2max(1,1),a1max(1,1),scalea1(1,1),CONSPA1)
      !A1x_half(i)=cal_pml_a(atype,Lx,L0,D2x_half(i),a1max(1,1),scalea1(1,1),CONSPA1)
      A2x_half(i)=cal_pml_a(a2type,Lx,L0, d1max(1,1),a2max(1,1),scalea2(1,1),CONSPA2)
    end do
  end if
  ! PML layer toward +x direction, starting from regular point where Tzz located
  if (abs_number(1,2)>0) then
    isx2=.true.
    x0=x(ni2-abs_number(1,2))
    L0=(x(ni2)+x(ni2+1))/2.0_SP-x0
    d1max(1,2)=cal_pml_dmax(L0,Vs(1,2),Rpp(1,2),b1max(1,2),CONSPD1,scaled1(1,2),a2max(1,2),d1type)
    d2max(1,2)=cal_pml_dmax(L0,Vs(1,2),Rpp(1,2),b2max(1,2),CONSPD2,scaled2(1,2),a1max(1,2),d2type)
    !  dmax(1,2)=0.0_SP
    do i=ni2-abs_number(1,2),ni2
      ! regular grid point
      Lx=x(i)-x0
      D1x_regu(i)=cal_pml_d(Lx,L0,d1max(1,2),CONSPD1)
      D2x_regu(i)=cal_pml_d(Lx,L0,d2max(1,2),CONSPD2)
      B1x_regu(i)=cal_pml_b(Lx,L0,b1max(1,2),CONSPB1)
      B2x_regu(i)=cal_pml_b(Lx,L0,b2max(1,2),CONSPB2)
      A1x_regu(i)=cal_pml_a(a1type,Lx,L0, d2max(1,2),a1max(1,2),scalea1(1,2),CONSPA1)
      !A1x_regu(i)=cal_pml_a(atype,Lx,L0,D2x_regu(i),a1max(1,2),scalea1(1,2),CONSPA1)
      A2x_regu(i)=cal_pml_a(a2type,Lx,L0, d1max(1,2),a2max(1,2),scalea2(1,2),CONSPA2)
      ! half grid point
      Lx=(x(i)+x(i+1))/2.0_SP-x0
      D1x_half(i)=cal_pml_d(Lx,L0,d1max(1,2),CONSPD1)
      D2x_half(i)=cal_pml_d(Lx,L0,d2max(1,2),CONSPD2)
      B1x_half(i)=cal_pml_b(Lx,L0,b1max(1,2),CONSPB1)
      B2x_half(i)=cal_pml_b(Lx,L0,b2max(1,2),CONSPB2)
      A1x_half(i)=cal_pml_a(a1type,Lx,L0, d2max(1,2),a1max(1,2),scalea1(1,2),CONSPA1)
      !A1x_half(i)=cal_pml_a(atype,Lx,L0,D2x_half(i),a1max(1,2),scalea1(1,2),CONSPA1)
      A2x_half(i)=cal_pml_a(a2type,Lx,L0, d1max(1,2),a2max(1,2),scalea2(1,2),CONSPA2)
    end do
  end if

  ! -------------------  z layer --------------------------
  ! PML layer toward -z direction, starting from half point where Vz located
  if (abs_number(2,1)>0) then
    isz1=.true.
    z0=(z(abs_number(2,1)+nk1+1)+z(abs_number(2,1)+nk1))/2.0_SP
    L0=z0-z(nk1)
    d1max(2,1)=cal_pml_dmax(L0,Vs(2,1),Rpp(2,1),b1max(2,1),CONSPD1,scaled1(2,1),a2max(2,1),d1type)
    d2max(2,1)=cal_pml_dmax(L0,Vs(2,1),Rpp(2,1),b2max(2,1),CONSPD2,scaled2(2,1),a1max(2,1),d2type)
    !  dmax(3,1)=0.0_SP
    do k=nk1,nk1+abs_number(2,1)
      ! regular grid point
      Lz=z0-z(k)
      D1z_regu(k)=cal_pml_d(Lz,L0,d1max(2,1),CONSPD1)
      D2z_regu(k)=cal_pml_d(Lz,L0,d2max(2,1),CONSPD2)
      B1z_regu(k)=cal_pml_b(Lz,L0,b1max(2,1),CONSPB1)
      B2z_regu(k)=cal_pml_b(Lz,L0,b2max(2,1),CONSPB2)
      A1z_regu(k)=cal_pml_a(a1type,Lz,L0, d2max(2,1),a1max(2,1),scalea1(2,1),CONSPA1)
      !A1z_regu(k)=cal_pml_a(atype,Lz,L0,D2z_regu(k),a1max(3,1),scalea1(3,1),CONSPA1)
      A2z_regu(k)=cal_pml_a(a2type,Lz,L0, d1max(2,1),a2max(2,1),scalea2(2,1),CONSPA2)
      ! half grid point
      Lz=z0-(z(k)+z(k+1))/2.0_SP
      D1z_half(k)=cal_pml_d(Lz,L0,d1max(2,1),CONSPD1)
      D2z_half(k)=cal_pml_d(Lz,L0,d2max(2,1),CONSPD2)
      B1z_half(k)=cal_pml_b(Lz,L0,b1max(2,1),CONSPB1)
      B2z_half(k)=cal_pml_b(Lz,L0,b2max(2,1),CONSPB2)
      A1z_half(k)=cal_pml_a(a1type,Lz,L0, d2max(2,1),a1max(2,1),scalea1(2,1),CONSPA1)
      !A1z_half(k)=cal_pml_a(atype,Lz,L0,D2z_half(k),a1max(3,1),scalea1(3,1),CONSPA1)
      A2z_half(k)=cal_pml_a(a2type,Lz,L0, d1max(2,1),a2max(2,1),scalea2(2,1),CONSPA2)
    end do
  end if

  ! PML layer toward +z direction, starting from regular point where Tzz located
  if (abs_number(2,2)>0) then
    isz2=.true.
    z0=z(nk2-abs_number(2,2))
    L0=(z(nk2)+z(nk2+1))/2.0_SP-z0
    d1max(2,2)=cal_pml_dmax(L0,Vs(2,2),Rpp(2,2),b1max(2,2),CONSPD1,scaled1(2,2),a2max(2,2),d1type)
    d2max(2,2)=cal_pml_dmax(L0,Vs(2,2),Rpp(2,2),b2max(2,2),CONSPD2,scaled2(2,2),a1max(2,2),d2type)
    !  dmax(3,2)=0.0_SP
    do k=nk2-abs_number(2,2),nk2
      ! regular grid point
      Lz=z(k)-z0
      D1z_regu(k)=cal_pml_d(Lz,L0,d1max(2,2),CONSPD1)
      D2z_regu(k)=cal_pml_d(Lz,L0,d2max(2,2),CONSPD2)
      B1z_regu(k)=cal_pml_b(Lz,L0,b1max(2,2),CONSPB1)
      B2z_regu(k)=cal_pml_b(Lz,L0,b2max(2,2),CONSPB2)
      A1z_regu(k)=cal_pml_a(a1type,Lz,L0, d2max(2,2),a1max(2,2),scalea1(2,2),CONSPA1)
      !A1z_regu(k)=cal_pml_a(atype,Lz,L0,D2z_regu(k),a1max(3,2),scalea1(3,2),CONSPA1)
      A2z_regu(k)=cal_pml_a(a2type,Lz,L0, d1max(2,2),a2max(2,2),scalea2(2,2),CONSPA2)
      ! half grkd point
      Lz=(z(k)+z(k+1))/2.0_SP-z0
      D1z_half(k)=cal_pml_d(Lz,L0,d1max(2,2),CONSPD1)
      D2z_half(k)=cal_pml_d(Lz,L0,d2max(2,2),CONSPD2)
      B1z_half(k)=cal_pml_b(Lz,L0,b1max(2,2),CONSPB1)
      B2z_half(k)=cal_pml_b(Lz,L0,b2max(2,2),CONSPB2)
      A1z_half(k)=cal_pml_a(a1type,Lz,L0, d2max(2,2),a1max(2,2),scalea1(2,2),CONSPA1)
      !A1z_half(k)=cal_pml_a(atype,Lz,L0,D2z_half(k),a1max(2,2),scalea1(2,2),CONSPA1)
      A2z_half(k)=cal_pml_a(a2type,Lz,L0, d1max(2,2),a2max(2,2),scalea2(2,2),CONSPA2)
    end do
  end if
  !== convert d_x to d_x/b_x ==
  !Dx_regu=Dx_regu/Bx_regu
  !Dz_regu=Dz_regu/Bz_regu
  !Dx_half=Dx_half/Bx_half
  !Dz_half=Dz_half/Bz_half

  if (isx1) then
    allocate( Vx11a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate( Vz11a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate(P11a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    Vx11a=0.0;  Vz11a=0.0
    P11a=0.0

    allocate( Vx12a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate( Vz12a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate(P12a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    Vx12a=0.0;  Vz12a=0.0
    P12a=0.0
  end if
  if (isx2) then
    allocate( Vx11b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate( Vz11b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate(P11b(ni2-abs_number(1,2):ni2,nk1:nk2))
    Vx11b=0.0;  Vz11b=0.0
    P11b=0.0

    allocate( Vx12b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate( Vz12b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate(P12b(ni2-abs_number(1,2):ni2,nk1:nk2))
    Vx12b=0.0;  Vz12b=0.0
    P12b=0.0
  end if

  if (isz1) then
    allocate( Vx31a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate( Vz31a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate(P31a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    Vx31a=0.0;  Vz31a=0.0
    P31a=0.0

    allocate( Vx32a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate( Vz32a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate(P32a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    Vx32a=0.0;  Vz32a=0.0
    P32a=0.0
  end if
  if (isz2) then
    allocate( Vx31b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate( Vz31b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate(P31b(ni1:ni2,nk2-abs_number(2,2):nk2))
    Vx31b=0.0;  Vz31b=0.0
    P31b=0.0

    allocate( Vx32b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate( Vz32b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate(P32b(ni1:ni2,nk2-abs_number(2,2):nk2))
    Vx32b=0.0;  Vz32b=0.0
    P32b=0.0
  end if

#ifdef DEBUG
  ! -- reflection coefficient --
  write(200+myid,*) "-------- reflection coefficient ----------"
  write(200+myid,"(a3,2es14.5)") 'x: ', Rpp(1,:)
  write(200+myid,"(a3,2es14.5)") 'z: ', Rpp(2,:)
  write(200+myid,"(a6,I5)") 'A1Tpye: ', a1type
  write(200+myid,"(a6,I5)") 'A2Tpye: ', a2type

  ! -- x damping --
  ! Ax at half and regular
  write(200+myid,*) "---- A1x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,A1x_half(i),A1x_regu(i)
  end do
  ! Bx at half and regular
  write(200+myid,*) "---- B1x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,B1x_half(i),B1x_regu(i)
  end do
  ! Dx at half and regular
  write(200+myid,*) "---- D1x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,D1x_half(i),D1x_regu(i)
  end do
  ! (Ax+Dx)*stept at half and regular
  write(200+myid,*) "---- (A1x(i)+D1x(i))*stept at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,(A1x_half(i)+D1x_half(i))*stept,(A1x_regu(i)+D1x_regu(i))*stept
  end do

  write(200+myid,*) "---- A2x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,A2x_half(i),A2x_regu(i)
  end do
  ! Bx at half and regular
  write(200+myid,*) "---- B2x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,B2x_half(i),B2x_regu(i)
  end do
  ! Dx at half and regular
  write(200+myid,*) "---- D2x(i) at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,D2x_half(i),D2x_regu(i)
  end do
  ! (Ax+Dx)*stept at half and regular
  write(200+myid,*) "---- (A2x(i)+D2x(i))*stept at half and regular points ----"
  do i=nx1,nx2
    write(200+myid,"(i6,2es14.5)") i,(A2x_half(i)+D2x_half(i))*stept,(A2x_regu(i)+D2x_regu(i))*stept
  end do

  ! -- z damping --
  ! Az at half and regular
  write(200+myid,*) "---- A1z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,A1z_half(k),A1z_regu(k)
  end do
  ! Bz at half and regular
  write(200+myid,*) "---- B1z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,B1z_half(k),B1z_regu(k)
  end do
  ! Dz at half and regular
  write(200+myid,*) "---- D1z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,D1z_half(k),D1z_regu(k)
  end do
  ! (Az+Dz)*stept at half and regular
  write(200+myid,*) "---- (A1z(k)+D1z(k))*stept at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,(A1z_half(k)+D1z_half(k))*stept,(A1z_regu(k)+D1z_regu(k))*stept
  end do

  write(200+myid,*) "---- A2z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,A2z_half(k),A2z_regu(k)
  end do
  ! Bz at half and regular
  write(200+myid,*) "---- B2z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,B2z_half(k),B2z_regu(k)
  end do
  ! Dz at half and regular
  write(200+myid,*) "---- D2z(k) at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,D2z_half(k),D2z_regu(k)
  end do
  ! (Az+Dz)*stept at half and regular
  write(200+myid,*) "---- (A2z(k)+D2z(k))*stept at half and regular points ----"
  do k=nz1,nz2
    write(200+myid,"(i6,2es14.5)") k,(A2z_half(k)+D2z_half(k))*stept,(A2z_regu(k)+D2z_regu(k))*stept
  end do

#endif

end subroutine abs_dpcfs_init

subroutine abs_dpcfs_destroy
  if (allocated(P11a)) deallocate(P11a)
  if (allocated(P31a)) deallocate(P31a)
  if (allocated( Vx11a)) deallocate( Vx11a)
  if (allocated( Vz11a)) deallocate( Vz11a)
  if (allocated( Vx31a)) deallocate( Vx31a)
  if (allocated( Vz31a)) deallocate( Vz31a)

  if (allocated(P12a)) deallocate(P11a)
  if (allocated(P32a)) deallocate(P31a)
  if (allocated( Vx12a)) deallocate( Vx11a)
  if (allocated( Vz12a)) deallocate( Vz11a)
  if (allocated( Vx32a)) deallocate( Vx31a)
  if (allocated( Vz32a)) deallocate( Vz31a)
end subroutine abs_dpcfs_destroy

!-----------------------------------------------------------------------------
subroutine get_para(Brx1,Brz1,  &
    Arx1,Arz1,  &
    Drx1,Drz1,  &
    Bhx1,Bhz1,  &
    Ahx1,Ahz1,  &
    Dhx1,Dhz1,  &
    Brx2,Brz2,  &
    Arx2,Arz2,  &
    Drx2,Drz2,  &
    Bhx2,Bhz2,  &
    Ahx2,Ahz2,  &
    Dhx2,Dhz2,  &
    Crx1,Crz1,  &
    Crx2,Crz2,  &
    Chx1,Chz1,  &
    Chx2,Chz2,  &
    i,   k)
  integer  ::    i,   k
  real(SP) :: Brx1,Brz1
  real(SP) :: Arx1,Arz1
  real(SP) :: Drx1,Drz1
  real(SP) :: Bhx1,Bhz1
  real(SP) :: Ahx1,Ahz1
  real(SP) :: Dhx1,Dhz1

  real(SP) :: Brx2,Brz2
  real(SP) :: Arx2,Arz2
  real(SP) :: Drx2,Drz2
  real(SP) :: Bhx2,Bhz2
  real(SP) :: Ahx2,Ahz2
  real(SP) :: Dhx2,Dhz2

  real(SP) :: Crx1,Crz1
  real(SP) :: Crx2,Crz2
  real(SP) :: Chx1,Chz1
  real(SP) :: Chx2,Chz2

  Brx1=B1x_regu(i); Bhx1=B1x_half(i)
  Arx1=A1x_regu(i); Ahx1=A1x_half(i)
  Drx1=D1x_regu(i); Dhx1=D1x_half(i)
  Brz1=B1z_regu(k); Bhz1=B1z_half(k)
  Arz1=A1z_regu(k); Ahz1=A1z_half(k)
  Drz1=D1z_regu(k); Dhz1=D1z_half(k)

  Brx2=B2x_regu(i); Bhx2=B2x_half(i)
  Arx2=A2x_regu(i); Ahx2=A2x_half(i)
  Drx2=D2x_regu(i); Dhx2=D2x_half(i)
  Brz2=B2z_regu(k); Bhz2=B2z_half(k)
  Arz2=A2z_regu(k); Ahz2=A2z_half(k)
  Drz2=D2z_regu(k); Dhz2=D2z_half(k)
  !Brx2=1.0_SP;   Bhx2=1.0_SP       
  !Arx2=0.0_SP;   Ahx2=0.0_SP
  !Drx2=0.1*Drx1; Dhx2=0.1*Dhx1
  !Brz2=1.0_SP;   Bhz2=1.0_SP
  !Arz2=0.0_SP;   Ahz2=0.0_SP
  !Drz2=0.1*Drz1; Dhz2=0.1*Dhz1

  !Arx1=cal_pml_a1(Drx2)
  !Ahx1=cal_pml_a1(Dhx2)
  !Arz1=cal_pml_a1(Drz2)
  !Ahz1=cal_pml_a1(Dhz2)

  Crx1=Drx1/Brx1*(1.0_SP-Drx2*Brx1/(Brx1*(Drx2+Arx2*Brx2)-Brx2*(Drx1+Arx1*Brx1)))
  Crx2=Drx2/Brx2*(1.0_SP-Drx1*Brx2/(Brx2*(Drx1+Arx1*Brx1)-Brx1*(Drx2+Arx2*Brx2)))
  Crz1=Drz1/Brz1*(1.0_SP-Drz2*Brz1/(Brz1*(Drz2+Arz2*Brz2)-Brz2*(Drz1+Arz1*Brz1)))
  Crz2=Drz2/Brz2*(1.0_SP-Drz1*Brz2/(Brz2*(Drz1+Arz1*Brz1)-Brz1*(Drz2+Arz2*Brz2)))

  Chx1=Dhx1/Bhx1*(1.0_SP-Dhx2*Bhx1/(Bhx1*(Dhx2+Ahx2*Bhx2)-Bhx2*(Dhx1+Ahx1*Bhx1)))
  Chx2=Dhx2/Bhx2*(1.0_SP-Dhx1*Bhx2/(Bhx2*(Dhx1+Ahx1*Bhx1)-Bhx1*(Dhx2+Ahx2*Bhx2)))
  Chz1=Dhz1/Bhz1*(1.0_SP-Dhz2*Bhz1/(Bhz1*(Dhz2+Ahz2*Bhz2)-Bhz2*(Dhz1+Ahz1*Bhz1)))
  Chz2=Dhz2/Bhz2*(1.0_SP-Dhz1*Bhz2/(Bhz2*(Dhz1+Ahz1*Bhz1)-Bhz1*(Dhz2+Ahz2*Bhz2)))
end subroutine get_para

subroutine abs_dpcfs_momentum
  integer :: n,i,k
  real(SP) :: DxP      
  real(SP) ::       DzP
  real(SP) :: P11,P31
  real(SP) :: P12,P32
  real(SP) :: rhox,rhoz

  real(SP) :: Brx1,Brz1
  real(SP) :: Arx1,Arz1
  real(SP) :: Drx1,Drz1
  real(SP) :: Bhx1,Bhz1
  real(SP) :: Ahx1,Ahz1
  real(SP) :: Dhx1,Dhz1

  real(SP) :: Brx2,Brz2
  real(SP) :: Arx2,Arz2
  real(SP) :: Drx2,Drz2
  real(SP) :: Bhx2,Bhz2
  real(SP) :: Ahx2,Ahz2
  real(SP) :: Dhx2,Dhz2

  real(SP) :: Crx1,Crz1
  real(SP) :: Crx2,Crz2
  real(SP) :: Chx1,Chz1
  real(SP) :: Chx2,Chz2
  if (isx1) then
    do k=nk1,nk2
    do i=ni1,ni1+abs_number(1,1)

      if(MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if

      DxP = m2d_FDxF(P,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct v_x
      P11 =  (2.0_SP * P11a(i,k) + stept*Chx1*DxP ) &
        /(2.0_SP + stept*(Ahx1+Dhx1/Bhx1))
      P12 =  (2.0_SP * P12a(i,k) + stept*Chx2*DxP ) &
        /(2.0_SP + stept*(Ahx2+Dhx2/Bhx2))
      Vx(i,k)=Vx(i,k)+stept/rhox*((1.0_SP/Bhx1/Bhx2-1)*DxP-P11/Bhx1/Bhx2-P12/Bhx1/Bhx2)
      P11a(i,k)=2.0_SP*P11-P11a(i,k)
      P12a(i,k)=2.0_SP*P12-P12a(i,k)

    end do
    end do
  end if

  if (isx2) then
    do k=nk1,nk2
    do i=ni2-abs_number(1,2),ni2

      if(MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if 

      DxP = m2d_FDxF(P,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct v_x
      P11 =  (2.0_SP * P11b(i,k) + stept*Chx1*DxP ) &
        /(2.0_SP + stept*(Ahx1+Dhx1/Bhx1))
      P12 =  (2.0_SP * P12b(i,k) + stept*Chx2*DxP ) &
        /(2.0_SP + stept*(Ahx2+Dhx2/Bhx2))
      Vx(i,k)=Vx(i,k)+stept/rhox*((1.0_SP/Bhx1/Bhx2-1)*DxP-P11/Bhx1/Bhx2-P12/Bhx1/Bhx2)
      P11b(i,k)=2.0_SP*P11-P11b(i,k)
      P12b(i,k)=2.0_SP*P12-P12b(i,k)

    end do
    end do
  end if

  if (isz1) then
    do k=nk1,nk1+abs_number(2,1)
    do i=ni1,ni2

      if (MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if

      DzP = m2d_FDzF(P,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct v_z
      P31 =  (2.0_SP * P31a(i,k) + stept*Chz1*DzP ) &
        /(2.0_SP + stept*(Ahz1+Dhz1/Bhz1))
      P32 =  (2.0_SP * P32a(i,k) + stept*Chz2*DzP ) &
        /(2.0_SP + stept*(Ahz2+Dhz2/Bhz2))
      Vz(i,k)=Vz(i,k)+stept/rhoz*((1.0_SP/Bhz1/Bhz2-1)*DzP -P31/Bhz1/Bhz2-P32/Bhz1/Bhz2)
      P31a(i,k)=2.0_SP*P31-P31a(i,k)
      P32a(i,k)=2.0_SP*P32-P32a(i,k)

    end do
    end do
  end if ! isz1

  if (isz2) then
    do k=nk2-abs_number(2,2),nk2
    do i=ni1,ni2

      if (MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if

      DzP = m2d_FDzF(P,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct v_z
      P31 =  (2.0_SP * P31b(i,k) + stept*Chz1*DzP ) &
        /(2.0_SP + stept*(Ahz1+Dhz1/Bhz1))
      P32 =  (2.0_SP * P32b(i,k) + stept*Chz2*DzP ) &
        /(2.0_SP + stept*(Ahz2+Dhz2/Bhz2))
      Vz(i,k)=Vz(i,k)+stept/rhoz*((1.0_SP/Bhz1/Bhz2-1)*DzP-P31/Bhz1/Bhz2-P32/Bhz1/Bhz2)
      P31b(i,k)=2.0_SP*P31-P31b(i,k)
      P32b(i,k)=2.0_SP*P32-P32b(i,k)

    end do
    end do
  end if ! isz2

end subroutine abs_dpcfs_momentum

subroutine abs_dpcfs_hook
  integer :: n,i,k
  real(SP) :: DxVx
  real(SP) :: DzVz
  real(SP) :: Vx11,Vx31,Vz11,Vz31
  real(SP) :: Vx12,Vx32,Vz12,Vz32
  real(SP) :: lam,miu,lam2mu

  real(SP) :: Brx1,Brz1
  real(SP) :: Arx1,Arz1
  real(SP) :: Drx1,Drz1
  real(SP) :: Bhx1,Bhz1
  real(SP) :: Ahx1,Ahz1
  real(SP) :: Dhx1,Dhz1

  real(SP) :: Brx2,Brz2
  real(SP) :: Arx2,Arz2
  real(SP) :: Drx2,Drz2
  real(SP) :: Bhx2,Bhz2
  real(SP) :: Ahx2,Ahz2
  real(SP) :: Dhx2,Dhz2

  real(SP) :: Crx1,Crz1
  real(SP) :: Crx2,Crz2
  real(SP) :: Chx1,Chz1
  real(SP) :: Chx2,Chz2
  if (isx1) then
    do k=nk1,nk2
    do i=ni1,ni1+abs_number(1,1)

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DxVx = m2d_FDxB(Vx,i,k)
      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct Txx,Tyy,Tzz
      Vx11 =  (2.0_SP * Vx11a(i,k) + stept*Crx1*DxVx ) &
        /(2.0_SP + stept*(Arx1+Drx1/Brx1))
      Vx12 =  (2.0_SP * Vx12a(i,k) + stept*Crx2*DxVx ) &
        /(2.0_SP + stept*(Arx2+Drx2/Brx2))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Brx1/Brx2-1)*DxVx-Vx11/Brx1/Brx2-Vx12/Brx1/Brx2)
      Vx11a(i,k)=2.0_SP*Vx11-Vx11a(i,k)
      Vx12a(i,k)=2.0_SP*Vx12-Vx12a(i,k)

      if (FreeOnTii .eqv. .true.) then
        if (freenode .and. k==nk2) then
          P(i,k)=P(i,k)-stept*lam**2/lam2mu*((1.0/Brx1/Brx2-1)*DxVx-Vx11/Brx1/Brx2-Vx12/Brx1/Brx2)
        end if
      end if

    end do
    end do
  end if ! isx1

  if (isx2) then
    do k=nk1,nk2
    do i=ni2-abs_number(1,2),ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;


      DxVx = m2d_FDxB(Vx,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct Txx,Tyy,Tzz
      Vx11 =  (2.0_SP * Vx11b(i,k) + stept*Crx1*DxVx ) &
        /(2.0_SP + stept*(Arx1+Drx1/Brx1))
      Vx12 =  (2.0_SP * Vx12b(i,k) + stept*Crx2*DxVx ) &
        /(2.0_SP + stept*(Arx2+Drx2/Brx2))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Brx1/Brx2-1)*DxVx-Vx11/Brx1/Brx2-Vx12/Brx1/Brx2)
      Vx11b(i,k)=2.0_SP*Vx11-Vx11b(i,k)
      Vx12b(i,k)=2.0_SP*Vx12-Vx12b(i,k)

      if (FreeOnTii .eqv. .true.) then
        if (freenode .and. k==nk2) then
          P(i,k)=P(i,k)-stept*lam**2/lam2mu*((1.0/Brx1/Brx2-1)*DxVx-Vx11/Brx1/Brx2-Vx12/Brx1/Brx2)
        end if
      end if

    end do
    end do
  end if ! isx2

  if (isz1) then
    do k=nk1,nk1+abs_number(2,1)
    do i=ni1,ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DzVz = m2d_FDzB(Vz,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct Txx,Tyy,Tzz
      Vz31 =  (2.0_SP * Vz31a(i,k) + stept*Crz1*DzVz ) &
        /(2.0_SP + stept*(Arz1+Drz1/Brz1))
      Vz32 =  (2.0_SP * Vz32a(i,k) + stept*Crz2*DzVz ) &
        /(2.0_SP + stept*(Arz2+Drz2/Brz2))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Brz1/Brz2-1)*DzVz-Vz31/Brz1/Brz2-Vz32/Brz1/Brz2)
      Vz31a(i,k)=2.0_SP*Vz31-Vz31a(i,k)
      Vz32a(i,k)=2.0_SP*Vz32-Vz32a(i,k)

    end do
    end do
  end if ! isz1

  if (isz2) then
    do k=nk2-abs_number(2,2),nk2
    do i=ni1,ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DzVz = m2d_FDzB(Vz,i,k)

      call get_para(Brx1,Brz1,  &
        Arx1,Arz1,  &
        Drx1,Drz1,  &
        Bhx1,Bhz1,  &
        Ahx1,Ahz1,  &
        Dhx1,Dhz1,  &
        Brx2,Brz2,  &
        Arx2,Arz2,  &
        Drx2,Drz2,  &
        Bhx2,Bhz2,  &
        Ahx2,Ahz2,  &
        Dhx2,Dhz2,  &
        Crx1,Crz1,  &
        Crx2,Crz2,  &
        Chx1,Chz1,  &
        Chx2,Chz2,  &
        i,   k)
      ! correct P
      Vz31 =  (2.0_SP * Vz31b(i,k) + stept*Crz1*DzVz ) &
        /(2.0_SP + stept*(Arz1+Drz1/Brz1))
      Vz32 =  (2.0_SP * Vz32b(i,k) + stept*Crz2*DzVz ) &
        /(2.0_SP + stept*(Arz2+Drz2/Brz2))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Brz1/Brz2-1)*DzVz-Vz31/Brz1/Brz2-Vz32/Brz1/Brz2)
      Vz31b(i,k)=2.0_SP*Vz31-Vz31b(i,k)
      Vz32b(i,k)=2.0_SP*Vz32-Vz32b(i,k)

    end do
    end do
  end if ! isz2

end subroutine abs_dpcfs_hook

!--------------------------------------------------------------------}

function cal_pml_dmax(L,Vp,Rpp,bmax,PD,scaled,amax,i) result(dmax)
  real(SP) :: L,Vp,Rpp,bmax,PD,scaled
  real(SP) :: dmax,amax
  integer  :: i
  if(i==1)then
    !dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)
    dmax=-Vp/2.0_SP/L*log(Rpp)*(PD+1.0_SP)*scaled
    !dmax=dmax*(CONSPD+CONSPB+1)/( (CONSPD+1.0_SP)*bmax+CONSPB )
  else
    dmax=amax*scaled
  endif
end function cal_pml_dmax

function cal_pml_d(x,L,dmax,PD) result(d)
  real(SP) :: x,L,dmax,PD
  real(SP) :: d
  !real(SP),parameter :: p=2.0_SP

  !dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)

  if (x<0) then
    d=0.0
  else
    d=dmax*(x/L)**PD
  end if
end function cal_pml_d

function cal_pml_amax(fc) result(amax)
  real(SP) :: fc
  real(SP) :: amax
  amax=PI*fc
end function cal_pml_amax

!function cal_pml_a(x,L,amax,PA) result(a)
!real(SP) :: x,L,amax,PA
!real(SP) :: a
!!amax=PI*fc
!if (x<0) then
!   a=0.0
!else
!   !a=amax*(L-x)/L
!   a=amax*(1.0_SP-(x/L)**PA)
!end if
!end function cal_pml_a

function cal_pml_a(i,x,L,d,amax,s,PA) result(a)
  real(SP) :: d,amax,x,L,s,PA
  integer  :: i
  real(SP) :: a
  if(i==1)then
    a=amax+d*s
  elseif(i==2)then
    a=amax
  elseif(i==3)then
    if (x<0) then
      a=0.0
    else
      a=amax*(1.0_SP-(x/L)**PA)+d*s
    end if
  elseif(i==4)then
    a=amax*(1.0_SP-(x/L)**PA)+d*(x/L)**2*s
  elseif(i==5)then
    a=amax+d*(x/L)**2*s
  endif
end function cal_pml_a

function cal_pml_b(x,L,bmax,PB) result(b)
  real(SP) :: x,L,bmax,PB
  real(SP) :: b
  !real(SP),parameter :: p=2.0_SP
  if (x<0) then
    b=1.0
  else
    b=1.0+(bmax-1.0)*(x/L)**PB
  end if
end function cal_pml_b

function cal_pml_R(n) result(r)
  integer :: n
  real(SP) :: r
  r = real(10**(-( log10(real(N,DP)) - 1.0_DP)/log10(2.0_DP) -3.0_DP),SP)
end function cal_pml_R
!---------------------- io related -----------------------------
subroutine abs_dpcfs_rest_export(pnm_rest)
  character (len=*) :: pnm_rest

end subroutine abs_dpcfs_rest_export

subroutine abs_dpcfs_rest_import(pnm_rest)
  character (len=*) :: pnm_rest


end subroutine abs_dpcfs_rest_import

end module abs_dpcfs_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
