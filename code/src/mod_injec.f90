

module injec_mod

! This module is used for seismic source injection
!
! Author: Wenzhong CAO  Email: caowz@mail.ustc.edu.cn
! Copyright (C) 2016 Wenzhong CAO

!*************************************************************************
!
! $Date: 2016-10-24 $
! $Revision: 001 $
! $LastChangedBy: caowz $
!
!*************************************************************************


































































use netcdf
use constants_mod
use string_mod, only : string_conf
use para_mod
use nompi_mod
use nfseis_mod
!use nfseis_mod, only : nfseis_diminfo,nfseis_except
!use io_mod, only : io_enum
use io_mod
use src_mod
use media_mod
use grid_mod
!use solver_mod 

implicit none
private

public ::                 &
  injec_src_read,         &
  injec_src_import,       &
  injec_src_hook,         &
  injec_src_momentum
!--     inner_wave_export
! src_injec_locate
!-- injec_src_locate

!------------------------------------------------------------------------! 
integer,public :: InjecSrc,injec_per_import,px1,pz1,px2,pz2

integer,allocatable,public :: src_box(:,:) 
integer,public :: injec_indx,injec_nt
real,allocatable,public :: injec_Txx(:,:),  &
  injec_Tzz(:,:),injec_Txz(:,:),            &
  injec_Vx(:,:), injec_Vz(:,:)
real,allocatable,public :: in_Txx(:,:),     &
  in_Tzz(:,:),in_Txz(:,:),                  &
  in_Vx(:,:),in_Vz(:,:)


character (len=SEIS_STRLEN),public :: fnm_src_injec,filenm
real,parameter,private :: MCC0=1.125, MCC1=0.0416666667

integer rest_tinv,run_from_rest

!------------------------------------------------------------------------! 
contains
!------------------------------------------------------------------------! 

!*************************************************************************
!*                       read source injection para                      *
!*************************************************************************

!!== read injection srouce parameters.
subroutine injec_src_read(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid,n,m
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
  call string_conf(fid,1,'InjecSrc',2,InjecSrc)
  write(*,*) "The source injection is ", InjecSrc

  if (InjecSrc/=1) return

  call alloc_injec_box

  do n=1,2
  do m=1,SEIS_GEO
    call string_conf(fid,1,trim(io_enum('point_',n)),m+1,src_box(m,n))
  end do
  end do
  call string_conf(fid,1,'injec_src_file'  ,2,fnm_src_injec)
  call string_conf(fid,1,'injec_per_import',2,injec_per_import)
  filenm=trim(pnm_src)//'/'//trim(fnm_src_injec)
  call nfseis_diminfo(filenm,'indx',injec_indx)
  call nfseis_diminfo(filenm,'nt',injec_nt)
  if (injec_per_import .eq. -1) then
    injec_per_import = injec_nt
  end if
  call alloc_injec_var

  px1=src_box(1,1)
  px2=src_box(1,2)
  pz1=src_box(2,1)
  pz2=src_box(2,2)

  call alloc_inner_var

  !-- write(*,*) "the injec parameters"
  write(*,*) "import ",injec_per_import," time steps per time"
  write(*,*) "point_001 is",src_box(:,1)
  write(*,*) "point_002 is",src_box(:,2)
  write(*,*) "fnm_src_injec file is     ",trim(fnm_src_injec)
  write(*,*) "filenm is     ",trim(filenm)
  write(*,*) "indx is ",injec_indx,"."," nt is ",injec_nt

  close(fid)
end subroutine injec_src_read

subroutine injec_src_hook(Txx,Tzz,Txz,ntime,stept)
  real(SP),dimension(:,:),intent(inout) :: Txx,Tzz,Txz
  integer,intent(in) :: ntime
  real(SP),intent(in) :: stept

  integer :: i,k,indxh

  !!== add injection srouce of outter area.
  indxh=0
  do k=pz1,pz2
  do i=px1,px2
    if ( k>(pz1+3) .and. k<(pz2-3) .and.  &
      i>(px1+3) .and. i<(px2-3) ) CYCLE
    indxh=indxh+1
    !!== inner area boundary
    !-- x direction.
    if (i<px1+2) then
      in_Txx(i,k) = Txx(i+2,k+2) - injec_Txx(indxh,ntime)
      in_Tzz(i,k) = Tzz(i+2,k+2) - injec_Tzz(indxh,ntime)
      in_Txz(i,k) = Txz(i+2,k+2) - injec_Txz(indxh,ntime)
      CYCLE
    end if
    if (i>px2-2) then
      in_Txx(i,k) = Txx(i+2,k+2) - injec_Txx(indxh,ntime)
      in_Tzz(i,k) = Tzz(i+2,k+2) - injec_Tzz(indxh,ntime)
      in_Txz(i,k) = Txz(i+2,k+2) - injec_Txz(indxh,ntime)
      CYCLE
    end if
    !-- z direction.
    if (k<pz1+2) then
      in_Txx(i,k) = Txx(i+2,k+2) - injec_Txx(indxh,ntime)
      in_Tzz(i,k) = Tzz(i+2,k+2) - injec_Tzz(indxh,ntime)
      in_Txz(i,k) = Txz(i+2,k+2) - injec_Txz(indxh,ntime)
      CYCLE
    end if
    if (k>pz2-2) then
      in_Txx(i,k) = Txx(i+2,k+2) - injec_Txx(indxh,ntime)
      in_Tzz(i,k) = Tzz(i+2,k+2) - injec_Tzz(indxh,ntime)
      in_Txz(i,k) = Txz(i+2,k+2) - injec_Txz(indxh,ntime)
      CYCLE
    end if
  end do ! px
  end do ! pz

  call cal_hook_inner(px1+2,px2-2,pz1+2,pz2-2)

  indxh=0
  do k=pz1,pz2
  do i=px1,px2
    if ( k>(pz1+3) .and. k<(pz2-3) .and.  &
      i>(px1+3) .and. i<(px2-3) ) CYCLE
    indxh=indxh+1
    !-- x direction.
    if ( i>px1+1 .and. i<px1+4 .and. &
      k>pz1+1 .and. k<pz2-1) then
      Txx(i+2,k+2)=injec_Txx(indxh,ntime) + in_Txx(i,k)
      Tzz(i+2,k+2)=injec_Tzz(indxh,ntime) + in_Tzz(i,k)
      Txz(i+2,k+2)=injec_Txz(indxh,ntime) + in_Txz(i,k)
      CYCLE
    end if
    if ( i>px2-4 .and. i<px2-1 .and. &
      k>pz1+1 .and. k<pz2-1) then
      Txx(i+2,k+2)=injec_Txx(indxh,ntime) + in_Txx(i,k)
      Tzz(i+2,k+2)=injec_Tzz(indxh,ntime) + in_Tzz(i,k)
      Txz(i+2,k+2)=injec_Txz(indxh,ntime) + in_Txz(i,k)
      CYCLE
    end if
    !-- z direction.
    if ( i>px1+1 .and. i<px2-1 .and. &
      k>pz1+1 .and. k<pz1+4) then
      Txx(i+2,k+2)=injec_Txx(indxh,ntime) + in_Txx(i,k)
      Tzz(i+2,k+2)=injec_Tzz(indxh,ntime) + in_Tzz(i,k)
      Txz(i+2,k+2)=injec_Txz(indxh,ntime) + in_Txz(i,k)
      CYCLE
    end if
    if ( i>px1+1 .and. i<px2-1 .and. &
      k>pz2-4 .and. k<pz2-1) then
      Txx(i+2,k+2)=injec_Txx(indxh,ntime) + in_Txx(i,k)
      Tzz(i+2,k+2)=injec_Tzz(indxh,ntime) + in_Tzz(i,k)
      Txz(i+2,k+2)=injec_Txz(indxh,ntime) + in_Txz(i,k)
      CYCLE
    end if
  end do ! px
  end do ! pz

  do k=pz1+4,pz2-4
  do i=px1+4,px2-4
    Txx(i+2,k+2) = in_Txx(i,k)
    Tzz(i+2,k+2) = in_Tzz(i,k)
    Txz(i+2,k+2) = in_Txz(i,k)
  end do
  end do

end subroutine injec_src_hook

subroutine injec_src_momentum(Vx,Vz,ntime,stept)
  real(SP),dimension(:,:),intent(inout) :: Vx,Vz
  integer,intent(in) :: ntime
  real(SP),intent(in) :: stept

  integer :: indxm,i,k

  !!== add injection srouce of outter area.
  indxm=0
  do k=pz1,pz2
  do i=px1,px2
    if ( k>(pz1+3) .and. k<(pz2-3) .and.  &
      i>(px1+3) .and. i<(px2-3) ) CYCLE
    indxm=indxm+1
    !-- x direction.
    if (i<px1+2) then
      in_Vx(i,k) = Vx(i+2,k+2) - injec_Vx(indxm,ntime)
      in_Vz(i,k) = Vz(i+2,k+2) - injec_Vz(indxm,ntime)
      CYCLE
    end if
    if (i>px2-2) then
      in_Vx(i,k) = Vx(i+2,k+2) - injec_Vx(indxm,ntime)
      in_Vz(i,k) = Vz(i+2,k+2) - injec_Vz(indxm,ntime)
      CYCLE
    end if
    !-- z direction.
    if (k<pz1+2) then
      in_Vx(i,k) = Vx(i+2,k+2) - injec_Vx(indxm,ntime)
      in_Vz(i,k) = Vz(i+2,k+2) - injec_Vz(indxm,ntime)
      CYCLE
    end if
    if (k>pz2-2) then
      in_Vx(i,k) = Vx(i+2,k+2) - injec_Vx(indxm,ntime)
      in_Vz(i,k) = Vz(i+2,k+2) - injec_Vz(indxm,ntime)
      CYCLE
    end if
  end do
  end do

  call cal_momentum_inner(px1+2,px2-2,pz1+2,pz2-2)

  indxm=0
  do k=pz1,pz2
  do i=px1,px2
    if ( k>(pz1+3) .and. k<(pz2-3) .and.  &
      i>(px1+3) .and. i<(px2-3) ) CYCLE
    indxm=indxm+1
    !-- x direction.
    if ( i>px1+1 .and. i<px1+4 .and. &
      k>pz1+1 .and. k<pz2-1) then
      Vx(i+2,k+2)=injec_Vx(indxm,ntime) + in_Vx(i,k)
      Vz(i+2,k+2)=injec_Vz(indxm,ntime) + in_Vz(i,k)
      CYCLE
    end if
    if ( i>px2-4 .and. i<px2-1 .and. &
      k>pz1+1 .and. k<pz2-1) then
      Vx(i+2,k+2)=injec_Vx(indxm,ntime) + in_Vx(i,k)
      Vz(i+2,k+2)=injec_Vz(indxm,ntime) + in_Vz(i,k)
      CYCLE
    end if
    !-- z direction.
    if ( i>px1+1 .and. i<px2-1 .and. &
      k>pz1+1 .and. k<pz1+4) then
      Vx(i+2,k+2)=injec_Vx(indxm,ntime) + in_Vx(i,k)
      Vz(i+2,k+2)=injec_Vz(indxm,ntime) + in_Vz(i,k)
      CYCLE
    end if
    if ( i>px1+1 .and. i<px2-1 .and. &
      k>pz2-4 .and. k<pz2-1) then
      Vx(i+2,k+2)=injec_Vx(indxm,ntime) + in_Vx(i,k)
      Vz(i+2,k+2)=injec_Vz(indxm,ntime) + in_Vz(i,k)
      CYCLE
    end if
  end do
  end do

  do k=pz1+4,pz2-4
  do i=px1+4,px2-4
    Vx(i+2,k+2) = in_Vx(i,k)
    Vz(i+2,k+2) = in_Vz(i,k)
  end do
  end do
end subroutine injec_src_momentum

!************************************************************************
!*                       update hook and momentum                       *
!************************************************************************
subroutine cal_hook_inner(I1,I2,K1,K2)
  integer,intent(in) :: I1,I2,K1,K2
  real(SP) :: lam,miu,lam2mu
  integer i,k
  real(SP) :: DxVx,DzVz
  !-- #ifdef CondFreeVAFDA
  !-- real(SP) :: DzVx
  !-- #endif
  integer :: k0

  !$OMP PARALLEL DEFAULT(shared)

  ! ========================  Tii ========================
  !$OMP DO PRIVATE(i,k,lam,miu,lam2mu)
  do k=K1,K2
  do i=I1,I2
    lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;
    DxVx=( MCC0*(in_Vx(i,k)-in_Vx(i-1,k)) -MCC1*(in_Vx(i+1,k)-in_Vx(i-2,k)) )/stepx
    DzVz=( MCC0*(in_Vz(i,k)-in_Vz(i,k-1)) -MCC1*(in_Vz(i,k+1)-in_Vz(i,k-2)) )/stepz
    in_Txx(i,k)=in_Txx(i,k) + stept*(      &
      lam2mu* DxVx + lam* DzVz )
    in_Tzz(i,k)=in_Tzz(i,k) + stept*(      &
      lam2mu* DzVz + lam* DxVx )
  end do
  end do
  !$OMP END DO

  ! ========================  Txz ========================
  !$OMP DO PRIVATE(i,k,miu)
  do k=K1,K2
  do i=I1,I2
    if (MediaPrecise .eqv. .true.) then
      miu=Uxz(i,k)
    else
      miu=4.0/(1.0/mu(i,k)+1.0/mu(i,k+1)+1.0/mu(i+1,k)+1.0/mu(i+1,k+1))
    end if
    in_Txz(i,k)=in_Txz(i,k) + stept*miu*( &
      ( MCC0*(in_Vz(i+1,k)-in_Vz(i,k)) -MCC1*(in_Vz(i+2,k)-in_Vz(i-1,k)) )/stepx             &
      +( MCC0*(in_Vx(i,k+1)-in_Vx(i,k)) -MCC1*(in_Vx(i,k+2)-in_Vx(i,k-1)) )/stepz )
  end do
  end do
  !$OMP END DO
end subroutine cal_hook_inner

subroutine cal_momentum_inner(I1,I2,K1,K2)
  integer,intent(in) :: I1,I2,K1,K2
  integer :: i,k
  real(SP) :: rhox,rhoz


  ! ========================  Vx ========================
  ! === interior domain upto nk2-2 ===
  !$OMP PARALLEL DEFAULT(shared)
  !$OMP DO PRIVATE(i,k,rhox)
  do k=K1,K2
  do i=I1,I2
    if(MediaPrecise .eqv. .true.) then
      rhox=Px(i,k)
    else
      rhox=0.5*(rho(i,k)+rho(i+1,k))
    end if
    in_Vx(i,k)=in_Vx(i,k) + stept*( &
      ( MCC0*(in_Txx(i+1,k)-in_Txx(i,k)) -MCC1*(in_Txx(i+2,k)-in_Txx(i-1,k)) )/stepx         &
      +( MCC0*(in_Txz(i,k)-in_Txz(i,k-1)) -MCC1*(in_Txz(i,k+1)-in_Txz(i,k-2)) )/stepz         &
      )/rhox
  end do
  end do
  !$OMP END DO

  ! ========================  Vz ========================
  !$OMP DO PRIVATE(i,j,k,rhoz)
  do k=K1,K2
  do i=I1,I2
    if(MediaPrecise .eqv. .true.) then
      rhoz=Pz(i,k)
    else
      rhoz=0.5*(rho(i,k)+rho(i,k+1))
    end if
    in_Vz(i,k)=in_Vz(i,k) + stept*( &
      ( MCC0*(in_Txz(i,k)-in_Txz(i-1,k)) -MCC1*(in_Txz(i+1,k)-in_Txz(i-2,k)) )/stepx         &
      +( MCC0*(in_Tzz(i,k+1)-in_Tzz(i,k)) -MCC1*(in_Tzz(i,k+2)-in_Tzz(i,k-1)) )/stepz         &
      )/rhoz
  end do
  end do
  !$OMP END DO

end subroutine cal_momentum_inner
!-- 
!-- 
!************************************************************************
!*                       netcdf format varget                           *
!************************************************************************

!!== import the injection sources
subroutine injec_src_import(ntime,nt)
  integer,intent(in) :: ntime,nt
  integer :: nimport
  integer,dimension(2) :: subs,subc,subt
  nimport = ntime/injec_per_import
  if (nt-ntime<injec_per_import) then
    subs=(/1,ntime+1/); subc=(/ injec_indx,nt-ntime+1 /); subt=(/1,1/)
  else
    subs=(/1,ntime+1/); subc=(/ injec_indx,injec_per_import /); subt=(/1,1/)
  end if

  print *, 'Injected source importing ',nimport
  print *, 'subs=',subs
  print *, 'subc=',subc
  !-- call injec_var_alloc
  call injec_varget(filenm,subs,subc,subt)
end subroutine injec_src_import

!!== netcdf format varget.
subroutine injec_varget(fnm,subs,subc,subt)
  character (len=*),intent(in) :: fnm
  integer,dimension(:),intent(in) :: subs,subc,subt
  integer :: ierr,ncid,Txxid,Tzzid,Txzid,   &
    Vxid,Vzid
  !-- open
  ierr = nf90_open(trim(fnm),NF90_NOWRITE,ncid)
  if (ierr /= nf90_noerr) & 
    call nfseis_except(ierr,'vargut_real1d:'//trim(fnm))

  !!== Txx
  ierr=nf90_inq_varid(ncid, 'Txx', Txxid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'var name in vargut_real1d:Txx')
  ierr=nf90_get_var(ncid,Txxid,injec_Txx,subs,subc,subt)
  if (ierr /= nf90_noerr) then
    print *, 'subs=',subs
    print *, 'subc=',subc
    print *, 'subt=',subt
    call nfseis_except(ierr,'vargut_real1d error:'//trim(fnm))
  end if

  !!== Tzz
  ierr=nf90_inq_varid(ncid, 'Tzz', Tzzid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'var name in vargut_real1d:Tzz')
  ierr=nf90_get_var(ncid,Tzzid,injec_Tzz,subs,subc,subt)
  if (ierr /= nf90_noerr) then
    print *, 'subs=',subs
    print *, 'subc=',subc
    print *, 'subt=',subt
    call nfseis_except(ierr,'vargut_real1d error:'//trim(fnm))
  end if

  !!== Txz
  ierr=nf90_inq_varid(ncid, 'Txz', Txzid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'var name in vargut_real1d:Txz')
  ierr=nf90_get_var(ncid,Txzid,injec_Txz,subs,subc,subt)
  if (ierr /= nf90_noerr) then
    print *, 'subs=',subs
    print *, 'subc=',subc
    print *, 'subt=',subt
    call nfseis_except(ierr,'vargut_real1d error:'//trim(fnm))
  end if

  !!== Vx
  ierr=nf90_inq_varid(ncid, 'Vx', Vxid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'var name in vargut_real1d:Vx')
  ierr=nf90_get_var(ncid,Vxid,injec_Vx,subs,subc,subt)
  if (ierr /= nf90_noerr) then
    print *, 'subs=',subs
    print *, 'subc=',subc
    print *, 'subt=',subt
    call nfseis_except(ierr,'vargut_real1d error:'//trim(fnm))
  end if

  !!== Vz
  ierr=nf90_inq_varid(ncid, 'Vz', Vzid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'var name in vargut_real1d:Vz')
  ierr=nf90_get_var(ncid,Vzid,injec_Vz,subs,subc,subt)
  if (ierr /= nf90_noerr) then
    print *, 'subs=',subs
    print *, 'subc=',subc
    print *, 'subt=',subt
    call nfseis_except(ierr,'vargut_real1d error:'//trim(fnm))
  end if

  !-- close
  ierr=nf90_close(ncid)
  if (ierr /= nf90_noerr) &
    call nfseis_except(ierr,'file close in vargut_real1d')
end subroutine injec_varget

!*************************************************************************
!*                       src alloc and dealloc                           *
!*************************************************************************
!!== allocate injec box, There are 2 point for 2D.
subroutine alloc_injec_box
  allocate(src_box(SEIS_GEO,2))
end subroutine alloc_injec_box

subroutine alloc_injec_var
  allocate(injec_Txx(injec_indx,injec_per_import)); injec_Txx=0.0
  allocate(injec_Tzz(injec_indx,injec_per_import)); injec_Tzz=0.0
  allocate(injec_Txz(injec_indx,injec_per_import)); injec_Txz=0.0
  allocate(injec_Vx (injec_indx,injec_per_import)); injec_Vx=0.0
  allocate(injec_Vz (injec_indx,injec_per_import)); injec_Vz=0.0
end subroutine alloc_injec_var

!!== allocate stress & vel.
subroutine alloc_inner_var
  allocate(in_Txx(px1:px2,pz1:pz2)); in_Txx=0.0
  allocate(in_Tzz(px1:px2,pz1:pz2)); in_Tzz=0.0
  allocate(in_Txz(px1:px2,pz1:pz2)); in_Txz=0.0
  allocate(in_Vx (px1:px2,pz1:pz2)); in_Vx=0.0
  allocate(in_Vz (px1:px2,pz1:pz2)); in_Vz=0.0
end subroutine alloc_inner_var

end module

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:

