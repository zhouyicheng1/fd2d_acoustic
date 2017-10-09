module src_mod

!----------------------------------------------------------------------------
! This module contains the varibales for 2D srouce.
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
use grid_mod
use media_mod
use nfseis_mod
implicit none
private

public ::           &
  src_fnm_init,     &
  src_alloc_moment, &
  src_alloc_force,  &
  src_alloc_resrc,  &
  src_destroy,      &
  src_import,       &
  src_fnm_get,      &
  src_resrc,        &
  src_stress,       &
  src_force
public ::        &
  fun_bell,      &
  fun_ricker,    &
  fun_gauss,     &
  src_stf,       &
  stf_name2id,   &
  stf_id2name

!----------------------------------------------------------------------!

integer,parameter ::      &
  SIG_SVF_BELL     =100,  &
  SIG_SVF_TRIANGLE =110,  &
  SIG_SVF_GAUSS    =120,  &
  SIG_STF_RICKER   =200,  &
  SIG_STF_BSHIFT   =210,  &
  SIG_STF_GAUSS    =220,  &
  SIG_STF_BELL     =230,  &
  SIG_STF_STEP     =240,  &
  SIG_STF_DELTA    =250

integer,public :: num_force,num_moment,num_resrc,ntwin_force,ntwin_moment
integer,allocatable,public :: force_indx(:,:),moment_indx(:,:),resrc_indx(:,:)
real(SP),allocatable,public ::                                  &
  force_shift(:,:), force_t0(:), frcstf_time(:),frcstf_freq(:), &
  moment_shift(:,:),moment_t0(:),momstf_time(:),momstf_freq(:)
real(SP),dimension(:,:),allocatable,public ::                   &
  resrc_Vx,resrc_Vz,                                            &
  MomTxx,MomTzz,MomTxy,MomTxz,MomTyz,                           &
  ForceX,ForceY,ForceZ
character (len=SEIS_STRLEN),public :: fnm_src_conf, pnm_src
integer,public :: momstf_id,frcstf_id,resrcstf_id
#ifdef VERBOSE
integer fid_out
#endif

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

!*************************************************************************
!*                ------ src alloc and dealloc ------                    *
!*************************************************************************
subroutine src_alloc_moment(npt,ntw)
  integer,intent(in) :: npt,ntw
  allocate(moment_indx(SEIS_GEO,npt)); moment_indx=0
  allocate(moment_shift(SEIS_GEO,npt)); moment_shift=0.0_SP
  allocate(moment_t0(npt)); moment_t0=0.0_SP
  allocate(MomTxx(ntw,npt)); MomTxx=0.0_SP
  allocate(MomTzz(ntw,npt)); MomTzz=0.0_SP
  allocate(MomTxz(ntw,npt)); MomTxz=0.0_SP
  allocate(MomTxy(ntw,npt)); MomTxy=0.0_SP
  allocate(MomTyz(ntw,npt)); MomTyz=0.0_SP
  allocate(momstf_time(ntw)); momstf_time=0.0_SP
  allocate(momstf_freq(ntw)); momstf_freq=0.0_SP
end subroutine src_alloc_moment
subroutine src_alloc_force(npt,ntw)
  integer,intent(in) :: npt,ntw
  allocate(force_indx(SEIS_GEO,npt)); force_indx=0
  allocate(force_shift(SEIS_GEO,npt)); force_shift=0.0_SP
  allocate(force_t0(npt)); force_t0=0.0_SP
  allocate(ForceX(ntw,npt)); ForceX=0.0_SP
  allocate(ForceZ(ntw,npt)); ForceZ=0.0_SP
  allocate(ForceY(ntw,npt)); ForceY=0.0_SP
  allocate(frcstf_time(ntw)); frcstf_time=0.0_SP
  allocate(frcstf_freq(ntw)); frcstf_freq=0.0_SP
end subroutine src_alloc_force
subroutine src_alloc_resrc(ntime,npt)
  integer,intent(in) :: npt,ntime
  allocate(resrc_indx(SEIS_GEO,npt)); resrc_indx=0
  allocate(resrc_Vx(ntw,npt)); resrc_Vx=0.0_SP
  allocate(resrc_Vz(ntw,npt)); resrc_Vz=0.0_SP
end subroutine src_alloc_resrc

subroutine src_destroy
!--- moment
  if (allocated(moment_indx)) deallocate(moment_indx)
  if (allocated(moment_shift)) deallocate(moment_shift)
  if (allocated(momstf_time)) deallocate(momstf_time)
  if (allocated(momstf_freq)) deallocate(momstf_freq)
  if (allocated(MomTxx)) deallocate(MomTxx)
  if (allocated(MomTzz)) deallocate(MomTzz)
  if (allocated(MomTxy)) deallocate(MomTxy)
  if (allocated(MomTxz)) deallocate(MomTxz)
  if (allocated(MomTyz)) deallocate(MomTyz)
!--- force
  if (allocated(force_indx)) deallocate(force_indx)
  if (allocated(force_shift)) deallocate(force_shift)
  if (allocated(frcstf_time)) deallocate(frcstf_time)
  if (allocated(frcstf_freq)) deallocate(frcstf_freq)
  if (allocated(ForceX)) deallocate(ForceX)
  if (allocated(ForceY)) deallocate(ForceY)
  if (allocated(ForceZ)) deallocate(ForceZ)
!--- resrc
  if (allocated(resrc_indx) deallocate(resrc_indx)
  if (allocated(resrc_Vx) deallocate(resrc_Vx)
  if (allocated(resrc_Vz) deallocate(resrc_Vz)
#ifdef VERBOSE
  close(fid_out)
#endif
end subroutine src_destroy

!*************************************************************************
!*                     ------ source io ------                           *
!*************************************************************************
subroutine src_fnm_init(filenm)
  character (len=*),intent(in) :: filenm
  integer fid
  fid=1001
  open(fid,file=trim(filenm),status="old")
  call string_conf(fid,1,'SOURCE_CONF',2,fnm_src_conf)
  call string_conf(fid,1,'SOURCE_ROOT',2,pnm_src)
  close(fid)
end subroutine src_fnm_init
function src_fnm_get() result(filenm)
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_src)//'/'//'source'//'.nc'
end function src_fnm_get

subroutine src_import
  integer n,i,k
  character (len=SEIS_STRLEN) :: filenm
  filenm=src_fnm_get()
  call nfseis_attget(filenm,'number_of_moment',num_moment)
  call nfseis_attget(filenm,'number_of_force',num_force)
  if (num_moment>0) then
    call nfseis_diminfo(filenm,'moment_time_window',ntwin_moment)
    call src_alloc_moment(num_moment,ntwin_moment)
    call nfseis_varget(filenm,'Mxx',MomTxx,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
    call nfseis_varget(filenm,'Mzz',MomTzz,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
    call nfseis_varget(filenm,'Mxz',MomTxz,(/1,1/),(/ntwin_moment,num_moment/),(/1,1/))
    call nfseis_varget(filenm,'moment_indx',moment_indx,(/1,1/),(/SEIS_GEO,num_moment/),(/1,1/))
    call nfseis_varget(filenm,'moment_shift',moment_shift,(/1,1/),(/SEIS_GEO,num_moment/),(/1,1/))
    call nfseis_varget(filenm,'moment_start_time',moment_t0,(/1/),(/num_moment/),(/1/))
    call nfseis_varget(filenm,'moment_rate_time',momstf_time,(/1/),(/ntwin_moment/),(/1/))
    call nfseis_varget(filenm,'moment_rate_freq',momstf_freq,(/1/),(/ntwin_moment/),(/1/))
    !att
    call nfseis_attget(filenm,'moment_rate_id',momstf_id)
    !correct
    !do n=1,num_moment
    !   moment_indx(:,n)=inn_ijk(moment_indx(:,n))
    !end do
  end if
  if (num_force>0) then
    call nfseis_diminfo(filenm,'force_time_window',ntwin_force)
    call src_alloc_force(num_force,ntwin_force)
    call nfseis_varget(filenm,'Fx',ForceX,(/1,1/),(/ntwin_force,num_force/),(/1,1/))
    call nfseis_varget(filenm,'Fz',ForceZ,(/1,1/),(/ntwin_force,num_force/),(/1,1/))

    call nfseis_varget(filenm,'force_indx',force_indx,(/1,1/),(/SEIS_GEO,num_force/),(/1,1/))
    call nfseis_varget(filenm,'force_shift',force_shift,(/1,1/),(/SEIS_GEO,num_force/),(/1,1/))
    call nfseis_varget(filenm,'force_start_time',force_t0,(/1/),(/num_force/),(/1/))
    call nfseis_varget(filenm,'force_stf_time',frcstf_time,(/1/),(/ntwin_force/),(/1/))
    call nfseis_varget(filenm,'force_stf_freq',frcstf_freq,(/1/),(/ntwin_force/),(/1/))
    !att
    call nfseis_attget(filenm,'force_stf_id',frcstf_id)
    !correct
    !do n=1,num_force
    !   force_indx(:,n)=inn_ijk(force_indx(:,n))
    !end do
  end if
  !resrc
  if (num_resrc>0) then
    call nfseis_diminfo(filenm,'time',ntt)
  !  call nfseis_diminfo(filenm,'resrc_point',num_resrc)
    call src_alloc_resrc(ntt,num_resrc)
    call nfseis_varget(filenm,'Vx',resrc_Vx,(/1,1/),(/ntt,num_resrc/),(/1,1/))
    call nfseis_varget(filenm,'Vz',resrc_Vz,(/1,1/),(/ntt,num_resrc/),(/1,1/))
    call nfseis_varget(filenm,'resrc_indx',resrc_indx,(/1,1/),(/SEIS_GEO,num_resrc/),(/1,1/))
  end if
end subroutine src_import

!*************************************************************************
!*                    ------ source coupling ------                      *
!*************************************************************************
subroutine src_stress(Txx,Tzz,Txz,ntime,stept)
 real(SP),dimension(:,:),intent(inout) :: Txx,Tzz,Txz
 integer, intent(in) :: ntime
 real(SP),intent(in) :: stept
 
 integer :: i,k,n,m,si,sk
 real(SP) :: t,rate, Mxx, Mxz, Mzz, V
 real(SP) :: dxx, dxz
 
!moment
 if (num_moment<=0) return
 t=cur_time
 do n=1,num_moment
    si=moment_indx(1,n);sk=moment_indx(2,n)
 do m=1,ntwin_moment
   rate=src_svf(t-moment_t0(n),momstf_time(m),momstf_freq(m),momstf_id)
   Mxx = MomTxx(m,n); Mzz=MomTzz(m,n); Mxz = MomTxz(m,n);
   i=si; k=sk; dxx=1.0_SP; dxz=1.0_SP;
   if (     i<=ni2 .and. i>=ni1    &
           .and. k<=nk2 .and. k>=nk1 ) then
       V=(x(i+1)-x(i-1))*(z(k+1)-z(k-1))/4.0
       V=rate*stept/V
       Txx(i,j,k)=Txx(i,j,k)-Mxx*V*dxx
       Tzz(i,j,k)=Tzz(i,j,k)-Mzz*V*dxx
       Txz(i,j,k)=Txz(i,j,k)-Mxz*V*dxz
    end if
end do
end do
end subroutine src_stress

subroutine src_force(Vx,Vz,ntime,stept)
real(SP),dimension(:,:),intent(inout) :: Vx,Vz
integer, intent(in) :: ntime
real(SP),intent(in) :: stept
integer :: i,k,n,m,si,sk
real(SP) :: t,disp,fx0,fz0,fx,fz,V
real(SP) :: dx,dz
real(SP) :: rhox,rhoz
!force
if ( num_force<=0)  return
t = cur_time
do n=1,num_force
   si=force_indx(1,n); sk=force_indx(2,n)
do m=1,ntwin_force
   disp=src_stf(t-force_t0(n),frcstf_time(m),frcstf_freq(m),frcstf_id)
   fx0=ForceX(m,n)*disp;  fz0=ForceZ(m,n)*disp
   i=si; k=sk; dx=1.0_SP;dz=1.0_SP
   if (      i<=ni2 .and. i>=ni1  &
       .and. k<=nk2 .and. k>=nk1 ) then
       rhox=Px(i,k); rhoz=Pz(i,k)
       V=(x(i+1)-x(i  ))*(z(k+1)-z(k-1))/4.0
       fx=dx*stept*fx0/rhox/V
       V=(x(i+1)-x(i-1))*(z(k+1)-z(k  ))/4.0
       fz=dz*stept*fz0/rhoz/V
       if (freenode .and. k==nk2) then
          fz=2.0_SP*fz
       end if
       Vx(i,k)=Vx(i,k)+fx
       Vz(i,k)=Vz(i,k)+fz
   end if
end do
end do
end subroutine src_force
subroutine src_resrc(Vx,Vz,Vax,Vaz,ntime,stept)
real(SP),dimension(:,:),intent(inout) :: Vx,Vz,Vax,Vaz
integer, intent(in) :: ntime
real(SP),intent(in) :: stept
real(SP) :: rx,rz
real(SP) :: rx0,rz0
real(SP) :: t,disp,V
integer :: i,k,n,m,si,sk
real(SP) :: rhox,rhoz
real(SP) :: dx,dz
!resrc
if ( num_resrc<=0)  return
t = cur_time
do n=1,num_resrc
   si=resrc_indx(1,n); sk=resrc_indx(2,n)
   rx =resrc_Vx(ntime,n); rz = resrc_Vz(ntime,n)
   i=si; j=sj; k=sk;
  if (      i<=ni2 .and. i>=ni1  &
      .and. k<=nk2 .and. k>=nk1 ) then

    Vx(i,k)=Vx(i,k)+rx
    Vz(i,k)=Vz(i,k)+rz
    Vax(i,k)=Vax(i,k)+rx
    Vaz(i,k)=Vaz(i,k)+rz
  end if
end do
end subroutine src_resrc


!#ifndef LITE
!subroutine src_surface_stress(TxSrc,TySrc,TzSrc,ntime)
!real,dimension(:,:) :: TxSrc,TySrc,TzSrc
!integer ntime
!
!integer i,j,k,n,m
!real t,rate,Mxz,Myz,Mzz
!real S
!
!if (flag_mech_type/=SIG_MECH_FORCE) return
!
!if (.not. freenode) return
!
!t=cur_time
!
!do n=1,npt_fault
!do m=1,ntwin
!   rate=cal_stf(t-RuptT(n)-twin(1,m), &
!                (twin(2,m)-twin(1,m)), stf_A(m) )
!   Mxz=MomTxx(m,n);Myz=MomTyy(m,n);Mzz=MomTzz(m,n)
!   i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)
!   if (      i<=ni2 .and. i>=ni1  &
!       .and. j<=nj2 .and. j>=nj1  &
!       .and. k==nk2 ) then
!       S=(x(i+1)-x(i-1))*(y(j+1)-y(j-1))/4.0
!       TzSrc(i,j)=Mzz*rate/S
!   end if
!   if (      i-1<=ni2 .and. i+1>=ni1  &
!       .and. j-1<=nj2 .and. j+1>=nj1  &
!       .and. k==nk2 ) then
!       S=(x(i  )-x(i-1))*(y(j+1)-y(j-1))/2.0
!       TxSrc(i-1,j)=Mxz*rate/S*( x(i+1)-x(i) )/( x(i+1)-x(i-1) )
!
!       S=(x(i+1)-x(i  ))*(y(j+1)-y(j-1))/2.0
!       TxSrc(i  ,j)=Mxz*rate/S*( x(i)-x(i-1) )/( x(i+1)-x(i-1) )
!
!       S=(x(i+1)-x(i-1))/2.0*(y(j  )-y(j-1))
!       TySrc(i,j-1)=Myz*rate/S*( y(j+1)-y(j) )/( y(j+1)-y(j-1) )
!
!       S=(x(i+1)-x(i-1))/2.0*(y(j+1)-y(j  ))
!       TySrc(i,j  )=Myz*rate/S*( y(j)-y(j-1) )/( y(j+1)-y(j-1) )
!   end if
!end do ! ntwin
!end do ! loop_of_npt_fault
!end subroutine src_surface_stress
!
!subroutine src_surface_velocity(VxSrc,VySrc,VzSrc,ntime)
!real,dimension(:,:) :: VxSrc,VySrc,VzSrc
!integer ntime
!integer i,j,k,n,m
!real t,rate,Mxz,Myz,Mzz
!real dx,dy,S
!
!t=cur_time
!
!if (.not. freenode) return
!
!do n=1,npt_fault
!do m=1,ntwin
!   rate=cal_svf(t-RuptT(n)-twin(1,m), &
!                (twin(2,m)-twin(1,m)), stf_A(m) )
!   Mxz=MomTxx(m,n);Myz=MomTyy(m,n);Mzz=MomTzz(m,n)
!   i=fault_indx(1,n); j=fault_indx(2,n); k=fault_indx(3,n)
!   if (      i<=ni2 .and. i>=ni1  &
!       .and. j<=nj2 .and. j>=nj1  &
!       .and. k==nk2 ) then
!       dx=(x(i+1)-x(i-1))/2.0
!       dy=(y(j+1)-y(j-1))/2.0
!       S=dx*dy
!       VzSrc(i,j)=Mzz*rate/S
!   end if
!   if (      i-1<=ni2 .and. i+1>=ni1  &
!       .and. j-1<=nj2 .and. j+1>=nj1  &
!       .and. k==nk2 ) then
!       S=(x(i  )-x(i-1))*(y(j+1)-y(j-1))/2.0
!       VxSrc(i-1,j)=Mxz*rate/S*( x(i+1)-x(i) )/( x(i+1)-x(i-1) )
!
!       S=(x(i+1)-x(i  ))*(y(j+1)-y(j-1))/2.0
!       VxSrc(i  ,j)=Mxz*rate/S*( x(i)-x(i-1) )/( x(i+1)-x(i-1) )
!
!       S=(x(i+1)-x(i-1))/2.0*(y(j  )-y(j-1))
!       VySrc(i,j-1)=Myz*rate/S*( y(j+1)-y(j) )/( y(j+1)-y(j-1) )
!
!       S=(x(i+1)-x(i-1))/2.0*(y(j+1)-y(j  ))
!       VySrc(i,j  )=Myz*rate/S*( y(j)-y(j-1) )/( y(j+1)-y(j-1) )
!   end if
!end do ! ntwin
!end do ! loop_of_npt_fault
!end subroutine src_surface_velocity
!#endif

!*************************************************************************
!*                --- source signal generator --                         *
!*************************************************************************
subroutine cal_norm_delt(delt,x0,z0,rx0,rz0)
  real(SP),dimension(-LenFD:LenFD,-LenFD:LenFD),intent(out) :: delt
  real(SP),intent(in) :: x0,z0,rx0,rz0
  real(SP) :: D1,D3
  integer i,k
  delt=0.0_SP
  do k=-LenFD,LenFD
  do i=-LenFD,LenFD
    D1=fun_gauss(i-x0,rx0,0.0_SP)
    D3=fun_gauss(k-z0,rz0,0.0_SP)
    delt(i,k)=D1*D3
  end do
  end do
  if (sum(delt)<SEIS_ZERO) call swmpi_except('cal_norm_delt is zero')
  delt=delt/sum(delt)
end subroutine cal_norm_delt
function src_svf(t,t0,f0,flag_stf_type) result(rate)
  real(SP),intent(in) :: t,t0,f0
  integer,intent(in) :: flag_stf_type
  real(SP) :: rate
  select case(flag_stf_type)
  case (SIG_SVF_GAUSS)
    rate=fun_gauss(t,f0,t0)
  case (SIG_STF_GAUSS)
    rate= fun_gauss_deriv(t,f0,t0)
  case (SIG_STF_RICKER)
    rate=fun_ricker_deriv(t,f0,t0)
  case (SIG_SVF_BELL)
    rate= fun_bell(t-t0,f0)
  case (SIG_STF_BELL)
    rate= fun_bell_deriv(t-t0,f0)
  case (SIG_SVF_TRIANGLE)
    rate= fun_triangle(t-t0,f0)
  case default
    print *, "Have you told me how to generate the SVF of ", flag_stf_type
    stop 1
  end select
  if (abs(rate)<SEIS_ZERO) rate=0.0
end function src_svf
function src_stf(t,t0,f0,flag_stf_type) result(disp)
  real(SP),intent(in) :: t,t0,f0
  integer,intent(in) :: flag_stf_type
  real(SP) :: disp
  select case(flag_stf_type)
  case (SIG_STF_GAUSS)
    disp=fun_gauss(t,f0,t0)
  case (SIG_STF_RICKER)
    disp=fun_ricker(t,f0,t0)
  case (SIG_STF_BELL)
    disp= fun_bell(t-t0,f0)
  case (SIG_SVF_BELL)
    disp= fun_bell_int(t-t0,f0)
  case (SIG_STF_STEP)
    disp=fun_step(t)
  case default
    print *, "Have you told me how to generate the STF of ", flag_stf_type
    stop 1
  end select
  !if (abs(disp)<SEIS_ZERO) disp=0.0
end function src_stf

! gauss
function fun_gauss(t,a,t0) result(f)
  real(SP),intent(in) :: t,t0,a
  real(SP) :: f
  if (abs(t0)>=SEIS_ZERO .and. (t<=0.0 .or. t>=2.0*t0) ) then
    f=0.0
  else
    f=exp(-(t-t0)**2/(a**2))/(sqrt(PI)*a)
  end if
end function fun_gauss
function fun_gauss_deriv(t,a,t0) result(f)
  real(SP),intent(in) :: t,t0,a
  real(SP) :: f
  if (abs(t0)>=SEIS_ZERO .and. (t<=0.0 .or. t>=2.0*t0) ) then
    f=0.0
  else
    f=exp(-(t-t0)**2/(a**2))/(sqrt(PI)*a)*(-2*(t-t0)/a**2)
  end if
end function fun_gauss_deriv
! ricker
function fun_ricker(t,fc,t0) result (v)
  real(SP),intent(in) :: t,t0,fc
  real(SP) :: u,f0,v
  if (t<=0.0) then
    v=0.0; return
  end if
  f0=sqrt(PI)/2.0
  u=(t-t0)*2.0*PI*fc
  v=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
  !v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker
function fun_ricker_deriv(t,fc,t0) result (v)
  real(SP),intent(in) :: t,t0,fc
  real(SP) :: u,f0,v
  if (t<=0.0) then
    v=0.0; return
  end if
  f0=sqrt(PI)/2.0
  u=(t-t0)*2.0*PI*fc
  !f=(u**2.0/4.0-0.5)*exp(-u**2.0/4.0)*f0
  v=u*(1.5-u**2/4.0)*exp(-u**2.0/4.0)*f0*PI*fc
end function fun_ricker_deriv
!bell
function fun_bell(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: v
  if (t>0.0 .and. t<riset) then
    v=(1.0 - cos(2*PI*t/riset))/riset
  else
    v=0.0
  end if
end function fun_bell
function fun_bell_deriv(t,riset) result(v)
  real,intent(in) :: t,riset
  real :: v
  if (t>0.0 .and. t<riset) then
    v=2.0*PI*sin(2*PI*t/riset)/riset
  else
    v=0.0
  end if
end function fun_bell_deriv
function fun_bell_int(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: v
  if (t<=0.0) then
    v=0.0
  elseif (t<riset) then
    v=t/riset - sin(2*PI*t/riset)/(2.0*pi)
  else
    v=1.0
  end if
end function fun_bell_int
!triangle
function fun_triangle(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: t0,v
  t0=riset/2.0;
  if (t>riset) then
    v=0.0
  elseif (t>t0) then
    v=2.0/t0-t/(t0**2)
  elseif (t>0.0) then
    v=t/(t0**2)
  else
    v=0.0
  end if
end function fun_triangle
function fun_triangle_int(t,riset) result(v)
  real(SP),intent(in) :: t,riset
  real(SP) :: t0,v
  t0=riset/2.0;
  if (t>riset) then
    v=1.0
  elseif (t>t0) then
    v=-0.5*t**2/(t0**2)+2.0*t/t0-1.0
  elseif (t>0.0) then
    v=0.5*t**2/(t0**2)
  else
    v=0.0
  end if
end function fun_triangle_int
!bshift
function fun_bshift(t,riset,t0) result(v)
  real(SP),intent(in) :: t,t0,riset
  real(SP) :: bshift,v
  bshift=riset/2.0;
  if (t>0.0 .and. t<riset) then
    v=0.5/t0/(cosh((t-bshift)/t0)**2.0)
  else
    v=0.0
  end if
end function fun_bshift
function fun_step(t) result(v)
  real(SP),intent(in) :: t
  real(SP) :: v
  if (t>0.0 ) then
    v=1.0
  else
    v=0.0
  end if
end function fun_step
function fun_delta(t,t0) result(v)
  real(SP),intent(in) :: t,t0
  real(SP) :: v
  if (t>=0.0 .and. t<t0) then
    v=1.0_SP/t0
  else
    v=0.0_SP
  end if
end function fun_delta

function stf_name2id(stfname) result(id)
  character (len=*),intent(in) :: stfname
  integer :: id
  select case (trim(stfname))
  case ('gauss')
    id=SIG_STF_GAUSS
  case ('gauss_int')
    id=SIG_SVF_GAUSS
  case ('ricker')
    id=SIG_STF_RICKER
  case ('bell')
    id=SIG_STF_BELL
  case ('bell_int')
    id=SIG_SVF_BELL
  case ('triangle_int')
    id=SIG_SVF_TRIANGLE
  case ('B(shift)')
    id=SIG_STF_BSHIFT
  case ('step')
    id=SIG_STF_STEP
  case default
    print *, "Have you told me how to generate the STF of ", trim(stfname)
    stop 1
  end select
end function stf_name2id
function stf_id2name(id) result(stfname)
  integer,intent(in) :: id
  character (len=SEIS_STRLEN) :: stfname
  select case (id)
  case (SIG_STF_GAUSS)
    stfname='gauss'
  case (SIG_SVF_GAUSS)
    stfname='gauss_int'
  case (SIG_STF_RICKER)
    stfname='ricker'
  case (SIG_STF_BELL)
    stfname='bell'
  case (SIG_SVF_BELL)
    stfname='bell_int'
  case (SIG_SVF_TRIANGLE)
    stfname='triangle_int'
  case (SIG_STF_BSHIFT)
    stfname='B(shift)'
  case (SIG_STF_STEP)
    stfname='step'
  case default
    print *, "Cann't find the stf name for ", id
    stop 1
  end select
end function stf_id2name

end module src_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:

