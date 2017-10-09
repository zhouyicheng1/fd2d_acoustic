module media_mod

!----------------------------------------------------------------------------
! This module contains the variables for 2D medium
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------

use constants_mod
use string_mod
use math_mod
use para_mod
use nompi_mod
use nfseis_mod

implicit none
private
public ::         &
  media_fnm_init, &
  media_fnm_get,  &
  media_destroy,  &
  media_alloc,    &
  media_import
!-- public ::                &
!--   media_alloc_lambdamu,  &
!--   media_alloc_rho,       &
!--   media_alloc_Qs,        &
!--   media_alloc_PUxyz
!-- public ::                &
!--   media_import_lambdamu, &
!--   media_import_rho,      &
!--   media_import_Qs,       &
!--   media_import_PUxyz

!-----------------------------------------------------------------------------

real(SP),dimension(:,:),allocatable,public :: rho,mu,lambda,Qs
real(SP),dimension(:,:),allocatable,public :: Px,Pz,Uxz
real(SP),public :: QsF0,QsINF
character (len=SEIS_STRLEN),public :: fnm_media_conf, pnm_media
integer :: ierr

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

!*************************************************************************
!*                            alloc and dealloc                          *
!*************************************************************************
subroutine media_alloc
  allocate( mu(nx1:nx2,nz1:nz2),stat=ierr);  mu=0.0
  allocate( lambda(nx1:nx2,nz1:nz2),stat=ierr);  lambda=0.0
  allocate( rho(nx1:nx2,nz1:nz2),stat=ierr);  rho=0.0
  if (MediaPrecise .eqv. .true.) then
    allocate( Px(nx1:nx2,nz1:nz2),stat=ierr);  Px=0.0
    allocate( Pz(nx1:nx2,nz1:nz2),stat=ierr);  Pz=0.0
    allocate( Uxz(nx1:nx2,nz1:nz2),stat=ierr);  Uxz=0.0
  end if
#ifdef WITHQS
  allocate( Qs(nx1:nx2,nz1:nz2),stat=ierr);  Qs=0.0
#endif
end subroutine media_alloc
subroutine media_destroy
  if (allocated(mu)) deallocate( mu )
  if (allocated(lambda)) deallocate( lambda )
  if (allocated(rho)) deallocate( rho )
  if (allocated(Qs)) deallocate( Qs )
end subroutine media_destroy

!*************************************************************************
!*                             media io                                  *
!*************************************************************************
subroutine media_fnm_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'MEDIA_CONF',2,fnm_media_conf)
  call string_conf(fid,1,'MEDIA_ROOT',2,pnm_media)
  close(fid)
end subroutine media_fnm_init

function media_fnm_get() result(filenm)
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_media)//'/'//'media'//'.nc'
end function media_fnm_get

subroutine media_import
  character (len=SEIS_STRLEN) :: filenm
  integer,dimension(SEIS_GEO) :: subs,subc,subt
  subs=(/ nx1,nz1 /); subc=(/ nx,nz /); subt=(/ 1,1 /)
  filenm=media_fnm_get()
  call nfseis_varget( filenm, 'mu', mu, subs,subc,subt)
  call nfseis_varget( filenm, 'lambda', lambda, subs,subc,subt)
  if(MediaPrecise .eqv. .true.) then
    call nfseis_varget( filenm, 'Px', Px, subs,subc,subt)
    call nfseis_varget( filenm, 'Pz', Pz, subs,subc,subt)
    call nfseis_varget( filenm, 'Uxz', Uxz, subs,subc,subt)
  else
    call nfseis_varget( filenm, 'rho', rho, subs,subc,subt)
  end if
#ifdef WITHQS
  call nfseis_varget( filenm, 'Qs', Qs, subs,subc,subt)
  call nfseis_attget( filenm, 'QsF0', QsF0)
  call nfseis_attget( filenm, 'QsINF', QsINF)
#endif
end subroutine media_import

end module media_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
