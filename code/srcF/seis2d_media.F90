program seis2d_media

!----------------------------------------------------------------------------
! This program generates media.
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
use math_mod
use para_mod
use nompi_mod
use nfseis_mod
use grid_mod
use media_mod

implicit none

integer :: NGRIDX,NGRIDZ

type STRUCT_INTERFACE
  logical :: yes
  character (len=SEIS_STRLEN) :: fnm
  character (len=SEIS_STRLEN) :: filetype
  integer :: ni,nk
  real(SP),dimension(:),pointer :: x
  real(SP),dimension(:,:),pointer :: z
  real(SP),dimension(:,:),pointer :: h
  real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
  real(SP),dimension(:,:),pointer :: Qs
  real(SP) :: QsF0,QsINF
end type STRUCT_INTERFACE

type STRUCT_LAYERED
  logical :: yes
  character (len=SEIS_STRLEN) :: fnm
  character (len=SEIS_STRLEN) :: filetype
  integer :: ni,nk
  real(SP),dimension(:),pointer :: x
  real(SP),dimension(:,:),pointer :: d
  real(SP),dimension(:,:),pointer :: h
  real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
  real(SP),dimension(:,:),pointer :: Qs
  real(SP) :: QsF0,QsINF
end type STRUCT_LAYERED

type STRUCT_COMPOSITE
  logical :: yes
  character (len=SEIS_STRLEN) :: fnm
  character (len=SEIS_STRLEN) :: filetype
  integer :: ni,nk
  real(SP),dimension(:),pointer :: x
  real(SP),dimension(:,:),pointer :: h
  real(SP),dimension(:,:,:),pointer :: Vp,Vs,Dp
  real(SP),dimension(:,:,:),pointer :: Qs
  real(SP),dimension(:),pointer :: Vp_poly_d,Vs_poly_d,Dp_poly_d,Qs_poly_d
  real(SP) :: QsF0,QsINF
end type STRUCT_COMPOSITE

type STRUCT_VOLUME
  logical :: yes
  character (len=SEIS_STRLEN) :: fnm
  character (len=SEIS_STRLEN) :: filetype
  integer :: imax,kmax
  integer :: ni,nk
  real(SP) :: rz0
  real(SP),dimension(:),pointer :: gx,grz
  real(SP),dimension(:),pointer :: x,rz
  real(SP),dimension(:,:),pointer :: Vp,Vs,Dp
  real(SP),dimension(:,:),pointer :: Qs
  real(SP) :: QsF0,QsINF
end type STRUCT_VOLUME

type (STRUCT_INTERFACE) :: LA
type (STRUCT_LAYERED) :: L2D
type (STRUCT_COMPOSITE) :: LC
type (STRUCT_VOLUME) :: BV,PV

integer :: i,k
integer,dimension(SEIS_GEO) :: subs,subc,subt
character (len=SEIS_STRLEN) :: filenm

real(SP) :: dtmax,dtmaxVp,dtmaxL
integer,dimension(SEIS_GEO) :: dtindx,dtnode

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)
call para_basic_init(fnm_conf)

call swmpi_init(fnm_conf)
call para_init(fnm_conf)

call swmpi_set_gindx(0,0)

call grid_fnm_init(fnm_conf)
call grid_coord_alloc

call media_fnm_init(fnm_conf)
call media_alloc

!-----------------------------------------------------------------------------

print *, 'init media ...'
call init_media(fnm_media_conf)

dtmax=1.0e10

!allocate(topogrid(nx)); topogrid=0.0
!allocate(topoz(nx)); topoz=0.0

!-----------------------------------------------------------------------------

print *, 'calculate effective media ...'

call grid_coord_import

if (BV%yes) call volume_read(BV)
if (PV%yes) call volume_read(PV)

call effmedia_eval

call media_extend

print *, '  export media...'
filenm=media_fnm_get()
call media_export(filenm)

call media_stept_check

print *, "Maximum allowed time step is", dtmax
!-- write(*,"(a,2i5)") "located on", dtindx
print *, " Vp and dL are:", dtmaxVp,dtmaxL
if (dtmax<stept) then
  print *, "Serious Error: stept>dtmax", stept,dtmax
  !stop 1
end if

call media_destroy
call grid_destroy

!-----------------------------------------------------------
contains
!-----------------------------------------------------------

subroutine init_media(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  character (len=SEIS_STRLEN) :: str
  integer :: fid

  fid=1001
  open(fid,file=trim(fnm_conf),status="old")

  call string_conf(fid,1,'half_sample_point',2,NGRIDX)
  call string_conf(fid,1,'half_sample_point',3,NGRIDZ)

  L2D%yes=.false.;
  LA%yes=.false.; LC%yes=.false.;
  BV%yes=.false.; PV%yes=.false.;

  ! background model
  call string_conf(fid,1,'background_type',2,str)
  select case (trim(str))
  case ('interface')
    LA%yes=.true.
    call string_conf(fid,1,'background_format',2,LA%filetype)
    call string_conf(fid,1,'background_filename',2,LA%fnm)
    call interface_read(LA)
  case ('layered')
    L2D%yes=.true.
    call string_conf(fid,1,'background_format',2,L2D%filetype)
    call string_conf(fid,1,'background_filename',2,L2D%fnm)
    call layered_read(L2D)
  case ('composite')
    LC%yes=.true.
    call string_conf(fid,1,'background_format',2,LC%filetype)
    call string_conf(fid,1,'background_filename',2,LC%fnm)
    call composite_read(LC)
  case ('volume')
    BV%yes=.true.
    call string_conf(fid,1,'background_format',2,BV%filetype)
    call string_conf(fid,1,'background_filename',2,BV%fnm)
    call volume_init(BV)
  end select

  ! perturbed model
  call string_conf(fid,1,'perturbed_type',2,str)
  select case (trim(str))
  case ('volume')
    PV%yes=.true.
    call string_conf(fid,1,'perturbed_format',2,PV%filetype)
    call string_conf(fid,1,'perturbed_filename',2,PV%fnm)
    call volume_init(PV)
  case ('none')
  end select

  close(fid)

  !if (W%yes .and. L%yes) then
  !   print *, "cann't deal with poly and layered media at same time currently"
  !   stop 1
  !end if
end subroutine init_media

! -----------------------------  interface structrue ----------------------------
subroutine interface_read(L)
  type(STRUCT_INTERFACE) :: L
  integer :: fid,imax,kmax,i,k
  real(DP) :: d2m,v2ms,d2kgm3
  character (len=SEIS_STRLEN) :: str

  if (trim(L%filetype)/='ascii') then
    call error_except('only ascii input accepted for interface type')
  end if

  fid=2001
  open(fid,file=trim(L%fnm),status="old")
  call string_conf(fid,1,'number_of_interface',2,kmax)
  if (kmax<1) then
    call error_except('at least 1 layer should exist for interface type')
  end if

  call string_conf(fid,1,'distance2meter',2,d2m)
  call string_conf(fid,1,'velocity2m/s',2,v2ms)
  call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
#ifdef WITHQS
  call string_conf(fid,1,'QsF0',2,L%QsF0)
  call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif
  !call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
  call string_conf(fid,1,'horizontal_sampling',2,imax)
  L%ni=imax; L%nk=kmax
  allocate(L%x(imax))
  allocate(L%z(imax,kmax))
  allocate(L%h(imax,kmax-1))
  allocate(L%Vp(2,kmax))
  allocate(L%Vs(2,kmax))
  allocate(L%Dp(2,kmax))
  allocate(L%Qs(2,kmax))

  call string_conf(fid,1,'<anchor_media>',1,str)
  do k=1,kmax
#ifdef WITHQS
    read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k),L%Qs(:,k)
#else
    read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k)
#endif
  end do

  !L%Vp(1,1)=L%Vp(2,1); L%Vs(1,1)=L%Vs(2,1)
  !L%Dp(1,1)=L%Dp(2,1); L%Qs(1,1)=L%Qs(2,1)
  L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%Dp=L%Dp*d2kgm3

  ! check
  do k=1,kmax
  do i=1,2
    if (k==1 .and. i==1) then
      if (L%Dp(i,k)<0 .or. L%Vp(i,k)<0 .or. L%Vs(i,k)<0 ) then
        call error_except('model parameter should >=0 for (1,1)')
      end if
    elseif ( L%Dp(i,k)<=0 .or. L%Vp(i,k)<=0 .or. L%Vs(i,k)<0 ) then
      print *, 'media parameter is negtive'
      print *, 'Vp=',L%Vp
      print *, 'Vs=',L%Vs
      print *, 'Dp=',L%Dp
      call error_except('interface type read failed')
    elseif ( L%Vp(i,k)**2<=2.0*L%Vs(i,k)**2 ) then
      print *, 'Vp**2 < 2Vs**2'
      print *, 'Vp=',L%Vp
      print *, 'Vs=',L%Vs
      call error_except('interface type read failed')
    end if
  end do
  end do

  call string_conf(fid,1,'<anchor_interface>',1,str)
  do i=1,imax     
    read(fid,*) L%x(i),( L%z(i,k), k=1,kmax )
  end do

  do k=1,kmax-1
    L%h(:,k)=L%z(:,k)-L%z(:,k+1)
  end do

  do k=1,kmax-1
    if (any(L%h(:,k)<0.0)) then
      print *, 'thickness of layer',k,' should be not less than 0'
      print *, minloc(L%h(:,k))
      call error_except("thickness btween two interfaces shouldn' be less than 0")
    end if
  end do
  L%x=L%x*d2m; L%z=L%z*d2m; L%h=L%h*d2m

  close(fid)

#ifdef WITHQS
  QsF0=L%QsF0; QsINF=L%QsINF
#endif
end subroutine interface_read

! -----------------------------  layered structrue ----------------------------
subroutine layered_read(L)
  type(STRUCT_LAYERED) :: L
  integer :: fid,imax,kmax,i,k
  real(DP) :: d2m,v2ms,d2kgm3
  character (len=SEIS_STRLEN) :: str,layer_type

  if (trim(L%filetype)/='ascii') then
    call error_except('only ascii input accepted for layered type')
  end if

  fid=2001
  open(fid,file=trim(L%fnm),status="old")
  call string_conf(fid,1,'number_of_layer',2,kmax)
  if (kmax<1) then
    call error_except('at least 1 layer should exist for layered type')
  end if

  call string_conf(fid,1,'distance2meter',2,d2m)
  call string_conf(fid,1,'velocity2m/s',2,v2ms)
  call string_conf(fid,1,'density2kg/m^3',2,d2kgm3)
#ifdef WITHQS
  call string_conf(fid,1,'QsF0',2,L%QsF0)
  call string_conf(fid,1,'QsINF',2,L%QsINF)
#endif
  !call string_conf(fid,1,'layer_sealevel',2,layer_sealevel)
  call string_conf(fid,1,'horizontal_sampling',2,imax)
  L%ni=imax; L%nk=kmax
  allocate(L%x(imax))
  allocate(L%d(imax,kmax))
  allocate(L%h(imax,kmax))
  allocate(L%Vp(2,kmax))
  allocate(L%Vs(2,kmax))
  allocate(L%Dp(2,kmax))
  allocate(L%Qs(2,kmax))

  call string_conf(fid,1,'<anchor_media>',1,str)
  do k=1,kmax
#ifdef WITHQS
    read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k),L%Qs(:,k)
#else
    read(fid,*) L%Vp(:,k),L%Vs(:,k),L%Dp(:,k)
#endif
  end do

  L%Vp=L%Vp*v2ms; L%Vs=L%Vs*v2ms; L%Dp=L%Dp*d2kgm3

  ! check
  if ( any(L%Dp<=0) .or. any(L%Vp<=0) .or. any(L%Vs<0) ) then
    print *, 'media parameter is negtive'
    print *, 'Vp=',L%Vp
    print *, 'Vs=',L%Vs
    print *, 'Dp=',L%Dp
    call error_except('layered type read failed')
  end if
  if ( any(L%Vp**2<=2.0*L%Vs**2) ) then
    print *, 'Vp**2 < 2Vs**2'
    print *, 'Vp=',L%Vp
    print *, 'Vs=',L%Vs
    call error_except('layered type read failed')
  end if

  call string_conf(fid,1,'layer_meaning',2,layer_type)
  call string_conf(fid,1,'<anchor_layer>',1,str)
  select case (trim(layer_type))
  case ('depth')
    do i=1,imax     
      read(fid,*) L%x(i),( L%d(i,k), k=1,kmax )
    end do
    L%h(:,1)=L%d(:,1)
    do k=2,kmax
      L%h(:,k)=L%d(:,k)-L%d(:,k-1)
    end do
  case ('thickness')
    do i=1,imax     
      read(fid,*) L%x(i),( L%h(i,k), k=1,kmax )
    end do
    L%d(:,1)=L%h(:,1)
    do k=2,kmax
      L%d(:,k)=L%d(:,k-1)+L%h(:,k)
    end do
  case default
    call error_except("layer_meaning can only take depth or thickness")
  end select

  do k=1,kmax-1
    if (any(L%h(:,k)<0.0)) then
      print *, minloc(L%h(:,k))
      call error_except('thickness of layer should not be less than 0')
    end if
  end do
  if (any(L%h(:,kmax)<SEIS_ZERO)) then
    print *, minloc(L%h(:,kmax))
    call error_except('thickness of lowest layer should larger than 0')
  end if

  L%x=L%x*d2m; L%d=L%d*d2m; L%h=L%h*d2m

  close(fid)

#ifdef WITHQS
  QsF0=L%QsF0; QsINF=L%QsINF
#endif
end subroutine layered_read

! -----------------------------  composite structrue ----------------------------
subroutine composite_read(L)
  type(STRUCT_COMPOSITE) :: L
  integer :: imax,kmax

  if (trim(L%filetype)/='nc') then
    call error_except('only nc input accepted for composite type')
  end if

  call nfseis_diminfo(L%fnm,'x',imax)
  call nfseis_diminfo(L%fnm,'layer',kmax)

  L%ni=imax; L%nk=kmax
  allocate(L%x(imax))
  allocate(L%h(imax,kmax))
  allocate(L%Vp(2,imax,kmax))
  allocate(L%Vs(2,imax,kmax))
  allocate(L%Dp(2,imax,kmax))
  allocate(L%Vp_poly_d(kmax))
  allocate(L%Vs_poly_d(kmax))
  allocate(L%Dp_poly_d(kmax))
#ifdef WITHQS
  allocate(L%Qs(2,imax,kmax))
  allocate(L%Qs_poly_d(kmax))
#endif

  call nfseis_varget(L%fnm,'x',L%x,(/1/),(/imax/),(/1/))
  call nfseis_varget(L%fnm,'thickness',L%h,(/1,1/),(/imax,kmax/),(/1,1/) )
  call nfseis_varget(L%fnm,'Vp',L%Vp,(/1,1,1/),(/2,imax,kmax/),(/1,1,1/) )
  call nfseis_varget(L%fnm,'Vs',L%Vs,(/1,1,1/),(/2,imax,kmax/),(/1,1,1/) )
  call nfseis_varget(L%fnm,'rho',L%Dp,(/1,1,1/),(/2,imax,kmax/),(/1,1,1/) )
  call nfseis_varget(L%fnm,'Vp_poly_d',L%Vp_poly_d,(/1/),(/kmax/),(/1/) )
  call nfseis_varget(L%fnm,'Vs_poly_d',L%Vs_poly_d,(/1/),(/kmax/),(/1/) )
  call nfseis_varget(L%fnm,'rho_poly_d',L%Dp_poly_d,(/1/),(/kmax/),(/1/) )
#ifdef WITHQS
  call nfseis_varget(L%fnm,'Qs',L%Qs,(/1,1,1/),(/2,imax,kmax/),(/1,1,1/) )
  call nfseis_attget(L%fnm,'QsF0',L%QsF0)
  call nfseis_attget(L%fnm,'QsINF',L%QsINF)
  QsF0=L%QsF0; QsINF=L%QsINF
#endif

  if (any(L%h(:,1:kmax-1)<0.0)) then
    call error_except('thickness of layer should not be less than 0')
  end if
  if (any(L%h(:,kmax)<SEIS_ZERO)) then
    print *, minloc(L%h(:,kmax))
    call error_except('thickness of lowest layer should larger than 0')
  end if
end subroutine composite_read

! -----------------------------  volume structrue ----------------------------
subroutine volume_init(P)
  type(STRUCT_VOLUME) :: P
  integer :: imax,kmax

  if (trim(P%filetype)/='nc') then
    call error_except('only nc input accepted for volume type')
  end if

  call nfseis_diminfo(P%fnm,'x',imax)
  call nfseis_diminfo(P%fnm,'depth',kmax)
  allocate(P%gx(imax)); P%gx=0.0
  allocate(P%grz(kmax)); P%grz=0.0
  call nfseis_varget(P%fnm,'x',P%gx,(/1/),(/imax/),(/1/))
  call nfseis_varget(P%fnm,'depth2sealevel',P%grz,(/1/),(/kmax/),(/1/))
  call nfseis_attget(P%fnm,'sealevel',P%rz0)
#ifdef WITHQS
  call nfseis_attget(P%fnm,'QsF0',P%QsF0)
  call nfseis_attget(P%fnm,'QsINF',P%QsINF)
  QsF0=P%QsF0; QsINF=P%QsINF
#endif
  P%imax=imax; P%kmax=kmax
  imax=min(nx*max(2*NGRIDX,1),imax)
  kmax=min(nz*max(2*NGRIDZ,1),kmax)
  call volume_alloc(P,imax,kmax)
end subroutine volume_init

subroutine volume_alloc(P,imax,kmax)
  type(STRUCT_VOLUME) :: P
  integer,intent(in) :: imax,kmax

  if (associated(P%Vp)) then
    if (size(P%Vp,1)<imax .or. size(P%Vp,2)<kmax) then
      deallocate(P%x);deallocate(P%rz);
      deallocate(P%Vp); deallocate(P%Vs); deallocate(P%Dp)
#ifdef WITHQS
      deallocate(P%Qs)
#endif
    end if
  end if

  P%ni=imax;P%nk=kmax
  if (.not. associated(P%Vp)) then
    allocate(P%x(imax)); P%x=0.0
    allocate(P%rz(kmax)); P%rz=0.0
    allocate(P%Vp(imax,kmax)); P%Vp=0.0
    allocate(P%Vs(imax,kmax)); P%Vs=0.0
    allocate(P%Dp(imax,kmax)); P%Dp=0.0
#ifdef WITHQS
    allocate(P%Qs(imax,kmax)); P%Qs=0.0
#endif
  end if
end subroutine volume_alloc

subroutine volume_read(P)
  type(STRUCT_VOLUME) :: P
  real(SP) :: xmin,xmax
  integer :: i1,i2,k1,k2,imax,kmax
  integer :: indx(1)
  xmin=minval(x); xmax=maxval(x)
  indx=maxloc(P%gx,P%gx<=max(xmin,P%gx(1)));i1=max(indx(1),1)
  indx=minloc(P%gx,P%gx>=min(xmax,P%gx(P%imax)));i2=indx(1)
  k1=1; k2=P%kmax
  imax=i2-i1+1; kmax=k2-k1+1
  call volume_alloc(P,imax,kmax)

  P%x(1:imax)=P%gx(i1:i2); P%rz(1:kmax)=P%grz(k1:k2);
  call nfseis_varget( P%fnm,'Vp',P%Vp(1:imax,1:kmax), &
    (/i1,k1/),(/imax,kmax/),(/1,1/) )
  call nfseis_varget( P%fnm,'Vs',P%Vs(1:imax,1:kmax), &
    (/i1,k1/),(/imax,kmax/),(/1,1/) )
  call nfseis_varget( P%fnm,'rho',P%Dp(1:imax,1:kmax), &
    (/i1,k1/),(/imax,kmax/),(/1,1/) )
#ifdef WITHQS
  call nfseis_varget( P%fnm,'Qs',P%Qs(1:imax,1:kmax), &
    (/i1,k1/),(/imax,kmax/),(/1,1/) )
#endif
end subroutine volume_read

! ---------------------------  ztopo structrue ----------------------------
!function ztopo_eval(x0,y0) result(ztopo)
!real(SP),intent(in) :: x0,y0
!real(SP) :: ztopo
!integer :: i0,j0
!integer :: p(SEIS_GEO-1)
!
!p=minloc( (topox-x0)**2+(topoy-y0)**2 )
!i0=p(1); j0=p(2)
!ztopo=topoz(i0,j0)
!!if (i0==1) i0=2
!!if (j0==1) j0=2
!!if (i0==NTPI) i0=NTPI-1
!!if (j0==NTPJ) j0=NTPJ-1
!!ztopo=interp_2d( topox(i0-1:i0+1,j0-1:j0+1),topoy(i0-1:i0+1,j0-1:j0+1), &
!!      topoz(i0-1:i0+1,j0-1:j0+1),3,3,x0,y0 )
!end function ztopo_eval

!---------------------------------

subroutine effmedia_eval

  real(SP),dimension(-NGRIDX:NGRIDX) :: xvec
  real(SP),dimension(-NGRIDZ:NGRIDZ) :: zvec
  real(SP) :: dx1,dx2,dz1,dz2
  real(SP) :: x0,y0,z0,d0,x1,x2,z1,z2

  real(SP),dimension(2) :: lam, miu
  real(SP),dimension(2) :: Vp,Vs,Dp,dVp,dVs,dDp
  real(SP) :: rho0,mu0,lam0
  logical :: flag_water
#ifdef WITHQS
  real(SP),dimension(2) :: Qm
  real(SP) :: Qs0
#endif

  integer :: nsamp
  integer :: i,k,mi,mk,n,m
  integer :: i1,i2,k1,k2,k0,n1,n2
  integer :: indx(1)

  !layered
  real(SP) :: h,d
  real(SP) :: L1,L2
  !verpoly
  integer :: wi,wk

  nsamp=max(2*NGRIDX,1)*max(2*NGRIDZ,1)
  if (MediaPrecise .eqv. .true.) then
    Px =0.0; Pz =0.0; Uxz=0.0
  end if

  ! one more for Uxz
  do k=nk1-LenFD+1,nk2+LenFD-1
    print *, ' k=',k-nk1+1, ' of',nk
    do i=ni1-LenFD+1,ni2+LenFD-1

      rho0=0.0;mu0=0.0;lam0=0.0
#ifdef WITHQS
      Qs0=0.0
#endif
      flag_water=.false.

      dx1=(x(i)-x(i-1))/2.0/(NGRIDX+1); dx2=(x(i+1)-x(i))/2.0/(NGRIDX+1)
      dz1=(z(k)-z(k-1))/2.0/(NGRIDZ+1); dz2=(z(k+1)-z(k))/2.0/(NGRIDZ+1)

      xvec(0)=x(i); zvec(0)=z(k)
      !xvec(-NGRIDX:-1)=(/-NGRIDX:-1/)*dx1+x(i);
      !xvec(1:NGRIDX)  =(/ 1:NGRIDX /)*dx2+x(i)
      !zvec(-NGRIDZ:-1)=(/-NGRIDZ:-1/)*dz1+z(k);
      !zvec( 1:NGRIDZ )=(/ 1:NGRIDZ /)*dz2+z(k)
      xvec(-NGRIDX:-1)=(/(n,n=-NGRIDX,-1,1)/)*dx1+x(i);
      xvec(1:NGRIDX)  =(/(n,n= 1,NGRIDX ,1)/)*dx2+x(i)
      zvec(-NGRIDZ:-1)=(/(n,n=-NGRIDZ,-1,1)/)*dz1+z(k);
      zvec( 1:NGRIDZ )=(/(n,n= 1,NGRIDZ ,1)/)*dz2+z(k)

      do mk=-NGRIDZ,NGRIDZ
        if (NGRIDZ/=0 .and. mk==0) cycle
        do mi=-NGRIDX,NGRIDX
          if (NGRIDX/=0 .and. mi==0) cycle

          x0=xvec(mi); z0=zvec(mk);

          ! -------------   background model ---------------------------
          ! interface
          if (LA%yes) then
            call indx_locate_1d(x0,LA%x,i1,i2,x1)

            do n=LA%nk,1,-1

              z1=interp_1d(LA%x(i1:i2),LA%z(i1:i2,n),2,x1)
              if (n==LA%nk) then
                h=1.0e3
              else
                h =interp_1d(LA%x(i1:i2),LA%h(i1:i2,n),2,x1)
              end if

              if (h<=SEIS_ZERO .and. n>1) cycle

              if (h<=SEIS_ZERO .and. n==1) then
                if (LA%Vp(1,n)>SEIS_ZERO) then
                  n1=1;k1=1;n2=1;k2=1
                end if
                Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
                Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
                Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
#ifdef WITHQS
                Qm(1)=LA%Qs(n1,k1); Qm(2)=LA%Qs(n2,k2)
#endif
                exit
              end if

              if (abs(z1-z0)<=dz2/5.0) then  !  .and. h>SEIS_ZERO
                k2=n; n2=2
                if (n==1) then
                  if (LA%Vp(1,n)>SEIS_ZERO) then
                    k1=n; n1=1
                  else
                    k1=n; n1=2
                  end if
                end if
                do m=n-1,1,-1
                  h=interp_1d(LA%x(i1:i2),LA%h(i1:i2,m),2,x1)
                  if (h>SEIS_ZERO) then
                    k1=m+1; n1=1
                    exit
                  elseif (m==1) then
                    if (LA%Vp(1,m)>SEIS_ZERO) then
                      k1=m; n1=1
                    else
                      k1=n; n1=2
                    end if
                  end if
                end do
                Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
                Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
                Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
#ifdef WITHQS
                Qm(1)=LA%Qs(n1,k1); Qm(2)=LA%Qs(n2,k2)
#endif
                exit
              end if

              if (z0<z1) then
                if (n==LA%nk) then
                  Vp(:)=LA%Vp(2,n)
                  Vs(:)=LA%Vs(2,n)
                  Dp(:)=LA%Dp(2,n)
#ifdef WITHQS
                  Qm(:)=LA%Qs(2,n)
#endif
                else
                  L2=(z1-z0)/h; L1=1.0-L2
                  Vp(:)=LA%Vp(2,n)*L1+LA%Vp(1,n+1)*L2
                  Vs(:)=LA%Vs(2,n)*L1+LA%Vs(1,n+1)*L2
                  Dp(:)=LA%Dp(2,n)*L1+LA%Dp(1,n+1)*L2
#ifdef WITHQS
                  Qm(:)=LA%Qs(2,n)*L1+LA%Qs(1,n+1)*L2
#endif
                end if
                exit
              end if

              k1=n; k2=n; n1=2; n2=2

              if (n==1) then
                if (LA%Vp(1,n)>SEIS_ZERO) then
                  n1=1;k1=1;n2=1;k2=1
                end if
                Vp(1)=LA%Vp(n1,k1); Vp(2)=LA%Vp(n2,k2)
                Vs(1)=LA%Vs(n1,k1); Vs(2)=LA%Vs(n2,k2)
                Dp(1)=LA%Dp(n1,k1); Dp(2)=LA%Dp(n2,k2)
#ifdef WITHQS
                Qm(1)=LA%Qs(n1,k1); Qm(2)=LA%Qs(n2,k2)
#endif
                exit
              end if

            end do
          end if
          ! layered
          if (L2D%yes) then
            d=0.0
            call indx_locate_1d(x0,L2D%x,i1,i2,x1)
            do n=1,L2D%nk
              h=interp_1d(L2D%x(i1:i2),L2D%h(i1:i2,n),2,x1)

              if (h<=SEIS_ZERO) cycle

              d=d+h
              if (abs(d-d0)<=dz2/5.0) then
                if (n==L2D%nk) call error_except('please increase the thickness of the lowest layer')
                n1=n;
                do m=n+1,L2D%nk
                  h=interp_1d(L2D%x(i1:i2),L2D%h(i1:i2,m),2,x1)
                  if (h>SEIS_ZERO) then
                    n2=m
                    exit
                  end if
                end do
                Vp(1)=L2D%Vp(2,n1); Vp(2)=L2D%Vp(1,n2)
                Vs(1)=L2D%Vs(2,n1); Vs(2)=L2D%Vs(1,n2)
                Dp(1)=L2D%Dp(2,n1); Dp(2)=L2D%Dp(1,n2)
#ifdef WITHQS
                Qm(1)=L2D%Qs(2,n1); Qm(2)=L2D%Qs(1,n2)
#endif
                exit
              elseif (d>d0) then
                L1=(d-d0)/h; L2=1.0-L1
                Vp(:)=L2D%Vp(1,n)*L1+L2D%Vp(2,n)*L2
                Vs(:)=L2D%Vs(1,n)*L1+L2D%Vs(2,n)*L2
                Dp(:)=L2D%Dp(1,n)*L1+L2D%Dp(2,n)*L2
#ifdef WITHQS
                Qm(:)=L2D%Qs(1,n)*L1+L2D%Qs(2,n)*L2
#endif
                exit
              elseif (n==L2D%nk) then
                call error_except('please increase the thickness of the lowest layer')
              end if
            end do
          end if
          ! composite
          if (LC%yes) then
            d=0.0
            call indx_locate_1d(x0,LC%x,i1,i2,x1)
            do n=1,LC%nk
              h=interp_1d(LC%x(i1:i2),LC%h(i1:i2,n),2,x1)

              if (h<=SEIS_ZERO) cycle

              d=d+h
              if (abs(d-d0)<=dz2/5.0) then
                if (n==LC%nk) call error_except('please increase the thickness of the lowest layer')
                n1=n
                do m=n+1,LC%nk
                  h=interp_1d(LC%x(i1:i2),LC%h(i1:i2,m),2,x1)
                  if (h>SEIS_ZERO) then
                    n2=m
                    exit
                  end if
                end do
                Vp(1)=interp_1d(LC%x(i1:i2),LC%Vp(2,i1:i2,n1),2,x1)
                Vs(1)=interp_1d(LC%x(i1:i2),LC%Vs(2,i1:i2,n1),2,x1)
                Dp(1)=interp_1d(LC%x(i1:i2),LC%Dp(2,i1:i2,n1),2,x1)
                Vp(2)=interp_1d(LC%x(i1:i2),LC%Vp(1,i1:i2,n2),2,x1)
                Vs(2)=interp_1d(LC%x(i1:i2),LC%Vs(1,i1:i2,n2),2,x1)
                Dp(2)=interp_1d(LC%x(i1:i2),LC%Dp(1,i1:i2,n2),2,x1)
#ifdef WITHQS
                Qm(1)=interp_1d(LC%x(i1:i2),LC%Qs(2,i1:i2,n1),2,x1)
                Qm(2)=interp_1d(LC%x(i1:i2),LC%Qs(1,i1:i2,n2),2,x1)
#endif
                exit
              elseif (d>d0) then
                Vp(1)=interp_1d(LC%x(i1:i2),LC%Vp(1,i1:i2,n),2,x1)
                Vs(1)=interp_1d(LC%x(i1:i2),LC%Vs(1,i1:i2,n),2,x1)
                Dp(1)=interp_1d(LC%x(i1:i2),LC%Dp(1,i1:i2,n),2,x1)
                Vp(2)=interp_1d(LC%x(i1:i2),LC%Vp(2,i1:i2,n),2,x1)
                Vs(2)=interp_1d(LC%x(i1:i2),LC%Vs(2,i1:i2,n),2,x1)
                Dp(2)=interp_1d(LC%x(i1:i2),LC%Dp(2,i1:i2,n),2,x1)
#ifdef WITHQS
                Qm(1)=interp_1d(LC%x(i1:i2),LC%Qs(1,i1:i2,n),2,x1)
                Qm(2)=interp_1d(LC%x(i1:i2),LC%Qs(2,i1:i2,n),2,x1)
#endif
                Vp(:)=Vp(1)+(Vp(2)-Vp(1))/h**LC%Vp_poly_d(n)*(d0-(d-h))**LC%Vp_poly_d(n)
                Vs(:)=Vs(1)+(Vs(2)-Vs(1))/h**LC%Vs_poly_d(n)*(d0-(d-h))**LC%Vs_poly_d(n)
                Dp(:)=Dp(1)+(Dp(2)-Dp(1))/h**LC%Dp_poly_d(n)*(d0-(d-h))**LC%Dp_poly_d(n)
#ifdef WITHQS
                Qm(:)=Qm(1)+(Qm(2)-Qm(1))/h**LC%Dp_poly_d(n)*(d0-(d-h))**LC%Dp_poly_d(n)
#endif
                exit
              elseif (n==LC%nk) then
                call error_except('please increase the thickness of the lowest layer')
              end if
            end do
          end if
          ! volume
          if (BV%yes) then
            call indx_locate_1d(x0,BV%x(1:BV%ni),i1,i2,x1)
            call indx_locate_1d(BV%rz0-z0,BV%rz(1:BV%nk),k1,k2,z1)
            Vp(:)=interp_2d(BV%x(i1:i2),BV%rz(k1:k2), &
              BV%Vp(i1:i2,k1:k2),2,2,x1,z1)
            Vs(:)=interp_2d(BV%x(i1:i2),BV%rz(k1:k2), &
              BV%Vs(i1:i2,k1:k2),2,2,x1,z1)
            Dp(:)=interp_2d(BV%x(i1:i2),BV%rz(k1:k2), &
              BV%Dp(i1:i2,k1:k2),2,2,x1,z1)
#ifdef WITHQS
            Qm(:)=interp_2d(BV%x(i1:i2),BV%rz(k1:k2), &
              BV%Qs(i1:i2,k1:k2),2,2,x1,z1)
#endif
          end if

          ! -------------   perturbed model ---------------------------
          ! volume
          if (PV%yes) then
            call indx_locate_1d(x0,PV%x(1:PV%ni),i1,i2,x1)
            call indx_locate_1d(PV%rz0-z0,PV%rz(1:PV%nk),k1,k2,z1)
            dVp=interp_2d(PV%x(i1:i2),PV%rz(k1:k2), &
              PV%Vp(i1:i2,k1:k2),2,2,x1,z1)
            dVs=interp_2d(PV%x(i1:i2),PV%rz(k1:k2), &
              PV%Vs(i1:i2,k1:k2),2,2,x1,z1)
            dDp=interp_2d(PV%x(i1:i2),PV%rz(k1:k2), &
              PV%Dp(i1:i2,k1:k2),2,2,x1,z1)
            Vp=Vp*(1.0+dVp)
            Vs=Vs*(1.0+dVs)
            Dp=Dp*(1.0+dDp)
          end if

          ! to lame parameters
          miu=Dp*Vs*Vs
          lam=Vp*Vp*Dp - 2.0*miu

          !accumulate
          rho0=rho0+0.5*(Dp(1)+Dp(2))
          if (miu(1)<=SEIS_ZERO .or. miu(2)<=SEIS_ZERO) then
            flag_water=.true.
          else
            mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
          end if
          !mu0 = mu0+0.5*(1.0/miu(1)+1.0/miu(2))
          lam0=lam0+0.5*(1.0/lam(1)+1.0/lam(2))
#ifdef WITHQS
          Qs0=Qs0+0.5*(Qm(1)+Qm(2))
#endif

          if (MediaPrecise .eqv. .true.) then
            wi=0; wk=0
            if (mi<0) wi=-1
            if (mk<0) wk=-1

            Px(i+wi,k)=Px(i+wi,k)+0.5*(Dp(1)+Dp(2))
            Pz(i,k+wk)=Pz(i,k+wk)+0.5*(Dp(1)+Dp(2))
            if (miu(1)<=SEIS_ZERO .or. miu(2)<=SEIS_ZERO) then
              flag_water=.true.
              Uxz(i+wi,k+wk)=Uxz(i+wi,k+wk)+0;
            else
              Uxz(i+wi,k+wk)=Uxz(i+wi,k+wk)+0.5*(1.0/miu(1)+1.0/miu(2))
            end if
          end if

        end do !mi
      end do !mk

      rho(i,k)=rho0/nsamp
      if (flag_water) then
        mu(i,k)=0.0
      else
        mu(i,k)=nsamp/mu0
      end if
      lambda(i,k)=nsamp/lam0
#ifdef WITHQS
      Qs(i,k)=Qs0/nsamp
#endif

    end do !i
  end do !k

  if (MediaPrecise .eqv. .true.) then
    do k=nk1,nk2
    do i=ni1,ni2
      Px (i,k)=Px (i,k)/nsamp
      Pz (i,k)=Pz (i,k)/nsamp
      if (Uxz(i,k)<=SEIS_ZERO) then
        Uxz(i,k)=0
      else
        Uxz(i,k)=nsamp/Uxz(i,k)
      end if
    end do
    end do
  end if
end subroutine effmedia_eval

function eval_poly(c,d,x) result(f)
  real(SP),dimension(:),intent(in) :: c,d
  real(SP),intent(in) :: x
  real(SP) :: f
  integer :: nmax,n
  nmax=size(d)
  !f=pp(1)*x0**3+pp(2)*x0**2+pp(3)*x0+pp(4)
  f=0.0
  do n=1,nmax
    f=f+c(n) * (x**d(n))
  end do
end function eval_poly

subroutine media_extend
  call extend_equal(rho)
  call extend_equal(mu)
  call extend_equal(lambda)
  if (MediaPrecise .eqv. .true.) then
    call extend_equal(Px)
    call extend_equal(Pz)
    call extend_equal(Uxz)
  end if
#ifdef WITHQS
  call extend_equal(Qs)
#endif
end subroutine media_extend

subroutine extend_equal(w)
  real(SP),dimension(nx1:nx2,nz1:nz2),intent(inout) :: w
  integer i,k,n
  ! x1, x2
  do k=nz1,nz2
  do n=1,LenFD
    w(ni1-n,k)=w(ni1,k)
    w(ni2+n,k)=w(ni2,k)
  end do
  end do
  ! z1, z2
  do i=nx1,nx2
  do n=1,LenFD
    w(i,nk1-n)=w(i,nk1)
    w(i,nk2+n)=w(i,nk2)
  end do
  end do
end subroutine extend_equal

function dist_point2line(x0,z0,p1,p2) result (L)
  real(SP),intent(in) :: x0,z0
  real(SP),dimension(SEIS_GEO),intent(in) :: p1,p2
  real(SP) :: L
  real(SP) :: A,B,C
  A=p2(2)-p1(2)
  B=P1(1)-p2(1)
  C=p2(1)*p1(2)-p1(1)*p2(2)
  L=abs((A*x0+B*z0+C)/sqrt(A**2+B**2))
end function dist_point2line

subroutine indx_locate_1d(x0,x,i1,i2,x1)
  real(SP),intent(in) :: x0
  real(SP),dimension(:),intent(in) :: x
  integer,intent(out) :: i1,i2
  real(SP),intent(out) :: x1
  integer :: i0
  !logical,intent(in) :: backward
  !integer :: indx(1)
  !indx=minloc(x,x>=x0)
  !i2=max(indx(1),2)
  !i1=i2-1

  i1=1; i2=size(x); x1=x0
  if (x0<=x(1)) then
    x1=x(1); i1=1; i2=i1+1
  elseif (x0>=x(i2)) then
    x1=x(i2); i1=i2-1
  else
    do while (i2-i1>1)
      i0=(i2-i1)/2+i1
      if (x(i0)==x0) then
        i1=i0;i2=i0+1; exit
      elseif (x(i0)>x0) then
        i2=i0;
      else
        i1=i0
      end if
    end do
  end if
end subroutine indx_locate_1d

subroutine dealloc_media
  if (L2D%yes) then
    deallocate(L2D%x)
    deallocate(L2D%h)
    deallocate(L2D%Vp)
    deallocate(L2D%Vs)
    deallocate(L2D%Dp)
  end if
end subroutine dealloc_media

subroutine media_export(filenm)
  character (len=*),intent(in) :: filenm
  integer,dimension(SEIS_GEO) ::subs,subc,subt
  call nfseis_grid2d_skel(filenm,nx,nz,"media generated by seis2d_media" )
  subs=(/ nx1,nz1 /); subc=(/ nx,nz /); subt=(/ 1,1 /)
  call nfseis_grid2d_attput(filenm,        &
    subs,subc,subt,                   &
    ni1,ni2,nk1,nk2,ni,nk, &
    nx1,nx2,nz1,nz2,nx,nz, &
    ngi1,ngi2,ngk1,ngk2,    &
    ngx1,ngx2,ngz1,ngz2,    &
    (/ngi1,ngi2,ngk1,ngk2/) )
#ifdef WITHQS
  call nfseis_attput(filenm,'QsF0',QsF0)
  call nfseis_attput(filenm,'QsINF',QsINF)
  call nfseis_grid2d_addvar(filenm,'Qs')
  call nfseis_varput(filenm,'Qs',Qs,subs,subc,subt)
#endif
  call nfseis_grid2d_addvar(filenm,'rho')
  call nfseis_grid2d_addvar(filenm,'mu')
  call nfseis_grid2d_addvar(filenm,'lambda')
  call nfseis_varput(filenm,'rho',rho,subs,subc,subt)
  call nfseis_varput(filenm,'mu',mu,subs,subc,subt)
  call nfseis_varput(filenm,'lambda',lambda,subs,subc,subt)
  if (MediaPrecise .eqv. .true.) then
    call nfseis_grid2d_addvar(filenm,'Px')
    call nfseis_grid2d_addvar(filenm,'Pz')
    call nfseis_grid2d_addvar(filenm,'Uxz')
    call nfseis_varput(filenm,'Px',Px,subs,subc,subt)
    call nfseis_varput(filenm,'Pz',Pz,subs,subc,subt)
    call nfseis_varput(filenm,'Uxz',Uxz,subs,subc,subt)
  end if
end subroutine media_export

subroutine media_stept_check
  real(SP) :: Vp,dtLe,dtlocal,L
  integer :: i,k
  integer :: ii,kk

  do k=nk1,nk2
  do i=ni1,ni2
    Vp=sqrt( (lambda(i,k)+2.0*mu(i,k))/rho(i,k) )

    dtLe = 1.e20;
    do kk=-1,1
    do ii=-1,1
      if(ii/=0 .and. kk/=0) then 
        L =  dist_point2line(x(i),z(k),  &
          (/ x(i-ii),z(k   ) /), &
          (/ x(i   ),z(k-kk) /))
        dtLe = min(dtLe, L)
      end if
    end do
    end do

    dtlocal=0.8487/Vp * dtLe 
    if (dtlocal<dtmax) then
      dtmax=dtlocal
      dtindx=(/ i,k /)
      dtmaxVp=Vp
      dtmaxL=dtLe
    end if
  end do
  end do
  write(*,"(a,1f12.10,a,2i5,a,2f12.10)") "   dtmax is", dtmax,          &
    " located on", dtindx," with Vp and dL", dtmaxVp,dtmaxL
end subroutine media_stept_check

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
  print *, trim(msg)
  stop 1
end subroutine error_except

end program seis2d_media

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
