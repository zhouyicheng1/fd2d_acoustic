module grid_mod

!----------------------------------------------------------------------------
! This module contains the variables for grid geometry
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
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
use string_mod
use para_mod
use nompi_mod
use nfseis_mod

implicit none

private
public ::             &
  grid_fnm_init,      &
  grid_coordfnm_get,  &
  grid_coord_alloc,   &
  grid_coord_import,  &
  grid_destroy
#ifdef NONUNIFORMGRID
public ::             &
  grid_fdcoef_alloc,  &
  grid_fdcoef_import
#endif

!-----------------------------------------------------------------------------
real(SP),public :: stepx,stepz

character (len=SEIS_STRLEN),public ::     &
  fnm_grid_conf,                          &
  pnm_grid
real(SP),allocatable,public :: x(:),z(:)

integer,public ::                         &
  npt_Vx, npt_Vz
integer,allocatable,public ::             &
  findx(:,:),                             &
  indxVx(:,:), indxVz(:,:)
real(SP),allocatable,public ::            &
  valVx(:),valVz(:)
integer,allocatable,public ::             &
  NyVx(:,:), NzVx(:),                     &
  NxVz(:,:), NyVz(:,:)

real(SP),dimension(:,:),allocatable,public :: CXF,CXB,CZF,CZB

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine grid_fnm_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer :: fid
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'GRID_CONF',2,fnm_grid_conf)
  call string_conf(fid,1,'GRID_ROOT',2,pnm_grid)
  close(fid)
end subroutine grid_fnm_init

subroutine grid_destroy
  if (allocated(findx)) deallocate(findx)
  if (allocated(x)) deallocate(x)
  if (allocated(z)) deallocate(z)
end subroutine grid_destroy

!*************************************************************************
!*                                  coord                                *
!*************************************************************************
subroutine grid_coord_alloc
  allocate( x(nx)); x=0.0_SP
  allocate( z(nz)); z=0.0_SP
end subroutine grid_coord_alloc

function grid_coordfnm_get() result(filenm)
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_grid)//'/'//'coord.nc'
end function grid_coordfnm_get

subroutine grid_coord_import
  character (len=SEIS_STRLEN) :: filenm
  filenm=grid_coordfnm_get()
  call nfseis_varget(filenm,'x',x,(/1/),(/nx/),(/1/))
  call nfseis_varget(filenm,'z',z,(/1/),(/nz/),(/1/))
#ifndef NONUNIFORMGRID
  call nfseis_attget(filenm,'stepx',stepx)
  call nfseis_attget(filenm,'stepz',stepz)
#endif
end subroutine grid_coord_import

#ifdef NONUNIFORMGRID
!*************************************************************************
!*                             nonuniform grid                           *
!*************************************************************************
subroutine grid_fdcoef_alloc
  allocate( CXF(-LenFDS:LenFDL,nx1:nx2)); CXF=0.0
  allocate( CXB(-LenFDL:LenFDS,nx1:nx2)); CXB=0.0
  allocate( CZF(-LenFDS:LenFDL,nz1:nz2)); CZF=0.0
  allocate( CZB(-LenFDL:LenFDS,nz1:nz2)); CZB=0.0
end subroutine grid_fdcoef_alloc
subroutine grid_fdcoef_import
  character (len=SEIS_STRLEN) :: filenm

  filenm=grid_metricfnm_get()

  call nfseis_varget(filenm,'CXF',CXF,(/1,1/),(/2*LenFD,nx/),(/1,1/))
  call nfseis_varget(filenm,'CXB',CXB,(/1,1/),(/2*LenFD,nx/),(/1,1/))
  call nfseis_varget(filenm,'CZF',CZF,(/1,1/),(/2*LenFD,nz/),(/1,1/))
  call nfseis_varget(filenm,'CZB',CZB,(/1,1/),(/2*LenFD,nz/),(/1,1/))
end subroutine grid_fdcoef_import
#endif

end module grid_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
