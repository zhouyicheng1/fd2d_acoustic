

module grid_mod

!----------------------------------------------------------------------------
! This module contains the variables for grid geometry
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------






































































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
  call nfseis_attget(filenm,'stepx',stepx)
  call nfseis_attget(filenm,'stepz',stepz)
end subroutine grid_coord_import


end module grid_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
