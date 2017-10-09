

module abs_cfs_mod

!----------------------------------------------------------------------------
! This module is used for absorbing outgoing waves based on
! unsplit-field CFS ADE-PML (auxiliary differential equation
! complex frequncy shifted PML)
! the complex stretching variable is sx = bx + dx/(ax + iw)
! This module implements PML on six faces after the main loop
! (the main loop includes the PML layers. Here uses PML
! auxiliary variable to correct the wavefield.
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------
































































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
public ::               &
  abs_cfs_init,         &
  abs_cfs_destroy,      &
  abs_cfs_hook,         &
  abs_cfs_momentum,     &
  abs_cfs_rest_export,  &
  abs_cfs_rest_import

!-----------------------------------------------------------------------------
real,parameter,private :: MCC0=1.125, MCC1=0.0416666667 

!integer,parameter :: SP=kind(1.0)
real(SP),parameter :: CONSPD=2.0_SP
real(SP),parameter :: CONSPB=2.0_SP
real(SP),parameter :: CONSPA=1.0_SP

real(SP),dimension(:,:),allocatable ::   &
  Vx1a,Vz1a,P1a, &
  Vx1b,Vz1b,P1b
real(SP),dimension(:,:),allocatable ::   &
  Vx3a,Vz3a,P3a, &
  Vx3b,Vz3b,P3b

logical :: isx1,isx2,isz1,isz2

integer,dimension(SEIS_GEO,2) :: abs_number
real(SP),dimension(:),allocatable :: Ax_regu,Az_regu
real(SP),dimension(:),allocatable :: Ax_half,Az_half
real(SP),dimension(:),allocatable :: Bx_regu,Bz_regu
real(SP),dimension(:),allocatable :: Bx_half,Bz_half
real(SP),dimension(:),allocatable :: Dx_regu,Dz_regu
real(SP),dimension(:),allocatable :: Dx_half,Dz_half

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine abs_cfs_init(fnm_conf)
  use nompi_mod, only : absnode,freenode
  character (len=*) :: fnm_conf
  integer fid,n,m,i,k,ierr,npt
  real(SP),dimension(SEIS_GEO,2) :: Vs,fc,bmax,Rpp,dmax,amax
  real(SP) :: x0,z0,L0,Lx,Lz

  fid=1001
  abs_number=0
  Vs=0.0_SP
  Rpp=0.0_SP
  dmax=0.0_SP
  amax=0.0_SP
  bmax=0.0_SP

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
      call string_conf(fid,1,'CFS_bmax',(n-1)*2+m+1,bmax(n,m))
      call string_conf(fid,1,'CFS_amax',(n-1)*2+m+1,amax(n,m))
      if (bmax(n,m)<1.0) then
        print *, "n,m,CFS_bmax(n,m)=",n,m,bmax(n,m)
        call swmpi_except("CFS_bmax should be large or equal to 1")
      end if
      Rpp(n,m)=cal_pml_R(abs_number(n,m))
    end if
  end do
  end do
  close(fid)

  !=== if layer >0, unset freenode ==
  if (abs_number(2,2)>0 .and. freenode) then
    freenode=.false.
  end if

  allocate(Ax_regu(nx1:nx2)); Ax_regu=0.0_SP
  allocate(Bx_regu(nx1:nx2)); Bx_regu=1.0_SP
  allocate(Dx_regu(nx1:nx2)); Dx_regu=0.0_SP
  allocate(Az_regu(nz1:nz2)); Az_regu=0.0_SP
  allocate(Bz_regu(nz1:nz2)); Bz_regu=1.0_SP
  allocate(Dz_regu(nz1:nz2)); Dz_regu=0.0_SP
  allocate(Ax_half(nx1:nx2)); Ax_half=0.0_SP
  allocate(Bx_half(nx1:nx2)); Bx_half=1.0_SP
  allocate(Dx_half(nx1:nx2)); Dx_half=0.0_SP
  allocate(Az_half(nz1:nz2)); Az_half=0.0_SP
  allocate(Bz_half(nz1:nz2)); Bz_half=1.0_SP
  allocate(Dz_half(nz1:nz2)); Dz_half=0.0_SP

  !== the PML starts from regular grid, where normal stress defined ==

  ! -------------------  x layer --------------------------
  ! PML layer toward -x direction, starting from half point where Vx located
  if (abs_number(1,1)>0) then
    isx1=.true.
    x0=(x(abs_number(1,1)+ni1+1)+x(abs_number(1,1)+ni1))/2.0_SP
    L0=x0-x(ni1)
    dmax(1,1)=cal_pml_dmax(L0,Vs(1,1),Rpp(1,1),bmax(1,1))
    do i=ni1,abs_number(1,1)+ni1
      ! regular grid point
      Lx=x0-x(i)
      Dx_regu(i)=cal_pml_d(Lx,L0,dmax(1,1))
      Ax_regu(i)=cal_pml_a(Lx,L0,amax(1,1))
      Bx_regu(i)=cal_pml_b(Lx,L0,bmax(1,1))
      ! half grid point
      Lx=x0-(x(i)+x(i+1))/2.0_SP
      Dx_half(i)=cal_pml_d(Lx,L0,dmax(1,1))
      Ax_half(i)=cal_pml_a(Lx,L0,amax(1,1))
      Bx_half(i)=cal_pml_b(Lx,L0,bmax(1,1))
    end do
  end if
  ! PML layer toward +x direction, starting from regular point where Tzz located
  if (abs_number(1,2)>0) then
    isx2=.true.
    x0=x(ni2-abs_number(1,2))
    L0=(x(ni2)+x(ni2+1))/2.0_SP-x0
    dmax(1,2)=cal_pml_dmax(L0,Vs(1,2),Rpp(1,2),bmax(1,2))
    do i=ni2-abs_number(1,2),ni2
      ! regular grid point
      Lx=x(i)-x0
      Dx_regu(i)=cal_pml_d(Lx,L0,dmax(1,2))
      Ax_regu(i)=cal_pml_a(Lx,L0,amax(1,2))
      Bx_regu(i)=cal_pml_b(Lx,L0,bmax(1,2))
      ! half grid point
      Lx=(x(i)+x(i+1))/2.0_SP-x0
      Dx_half(i)=cal_pml_d(Lx,L0,dmax(1,2))
      Ax_half(i)=cal_pml_a(Lx,L0,amax(1,2))
      Bx_half(i)=cal_pml_b(Lx,L0,bmax(1,2))
    end do
  end if

  ! -------------------  z layer --------------------------
  ! PML layer toward -z direction, starting from half point where Vz located
  if (abs_number(2,1)>0) then
    isz1=.true.
    z0=(z(abs_number(2,1)+nk1+1)+z(abs_number(2,1)+nk1))/2.0_SP
    L0=z0-z(nk1)
    dmax(2,1)=cal_pml_dmax(L0,Vs(2,1),Rpp(2,1),bmax(2,1))
    do k=nk1,nk1+abs_number(2,1)
      ! regular grid point
      Lz=z0-z(k)
      Dz_regu(k)=cal_pml_d(Lz,L0,dmax(2,1))
      Az_regu(k)=cal_pml_a(Lz,L0,amax(2,1))
      Bz_regu(k)=cal_pml_b(Lz,L0,bmax(2,1))
      ! half grid point
      Lz=z0-(z(k)+z(k+1))/2.0_SP
      Dz_half(k)=cal_pml_d(Lz,L0,dmax(2,1))
      Az_half(k)=cal_pml_a(Lz,L0,amax(2,1))
      Bz_half(k)=cal_pml_b(Lz,L0,bmax(2,1))
    end do
  end if

  ! PML layer toward +z direction, starting from regular point where Tzz located
  if (abs_number(2,2)>0) then
    isz2=.true.
    z0=z(nk2-abs_number(2,2))
    L0=(z(nk2)+z(nk2+1))/2.0_SP-z0
    dmax(2,2)=cal_pml_dmax(L0,Vs(2,2),Rpp(2,2),bmax(2,2))
    do k=nk2-abs_number(2,2),nk2
      ! regular grid point
      Lz=z(k)-z0
      Dz_regu(k)=cal_pml_d(Lz,L0,dmax(2,2))
      Az_regu(k)=cal_pml_a(Lz,L0,amax(2,2))
      Bz_regu(k)=cal_pml_b(Lz,L0,bmax(2,2))
      ! half grkd point
      Lz=(z(k)+z(k+1))/2.0_SP-z0
      Dz_half(k)=cal_pml_d(Lz,L0,dmax(2,2))
      Az_half(k)=cal_pml_a(Lz,L0,amax(2,2))
      Bz_half(k)=cal_pml_b(Lz,L0,bmax(2,2))
    end do
  end if

  !== convert d_x to d_x/b_x ==
  Dx_regu=Dx_regu/Bx_regu
  Dz_regu=Dz_regu/Bz_regu
  Dx_half=Dx_half/Bx_half
  Dz_half=Dz_half/Bz_half

  if (isx1) then
    allocate( Vx1a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate( Vz1a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    allocate(P1a(ni1:abs_number(1,1)+ni1,nk1:nk2))
    Vx1a=0.0;  Vz1a=0.0
    P1a=0.0
  end if
  if (isx2) then
    allocate( Vx1b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate( Vz1b(ni2-abs_number(1,2):ni2,nk1:nk2))
    allocate(P1b(ni2-abs_number(1,2):ni2,nk1:nk2))
    Vx1b=0.0;  Vz1b=0.0
    P1b=0.0
  end if

  if (isz1) then
    allocate( Vx3a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate( Vz3a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    allocate(P3a(ni1:ni2,nk1:abs_number(2,1)+nk1))
    Vx3a=0.0;  Vz3a=0.0
    P3a=0.0
  end if
  if (isz2) then
    allocate( Vx3b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate( Vz3b(ni1:ni2,nk2-abs_number(2,2):nk2))
    allocate(P3b(ni1:ni2,nk2-abs_number(2,2):nk2))
    Vx3b=0.0;  Vz3b=0.0
    P3b=0.0
  end if


end subroutine abs_cfs_init

subroutine abs_cfs_destroy
  if (allocated(P1a)) deallocate(P1a)
  if (allocated(P3a)) deallocate(P3a)
  if (allocated( Vx1a)) deallocate( Vx1a)
  if (allocated( Vz1a)) deallocate( Vz1a)
  if (allocated( Vx3a)) deallocate( Vx3a)
  if (allocated( Vz3a)) deallocate( Vz3a)
end subroutine abs_cfs_destroy

!-----------------------------------------------------------------------------

subroutine abs_cfs_momentum
  integer :: n,i,k
  real(SP) :: DxP,DzP
  real(SP) :: P1,P3
  real(SP) :: rhox,rhoz

  if (isx1) then
    do k=nk1,nk2
    do i=ni1,ni1+abs_number(1,1)

      if (MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if

      DxP = ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx

      ! correct v_x
      P1 =  (2.0_SP * P1a(i,k) + stept*Dx_half(i)*DxP ) &
        /(2.0_SP + stept*(Ax_half(i)+Dx_half(i)))
      Vx(i,k)=Vx(i,k)+stept/rhox*((1.0_SP/Bx_half(i)-1)*DxP-P1/Bx_half(i))
      P1a(i,k)=2.0_SP*P1-P1a(i,k)

    end do
    end do
  end if

  if (isx2) then
    do k=nk1,nk2
    do i=ni2-abs_number(1,2),ni2

      if (MediaPrecise .eqv. .true.) then
        rhox=Px(i,k)
        rhoz=Pz(i,k)
      else
        rhox=0.5_sp*(rho(i,k)+rho(i+1,k))
        rhoz=0.5_sp*(rho(i,k)+rho(i,k+1))
      end if

      DxP = ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx
      !DxTxz = ( MCC0*(Txz(i,k)-Txz(i-1,k)) -MCC1*(Txz(i+1,k)-Txz(i-2,k)) )/stepx

      ! correct v_x
      P1 =  (2.0_SP * P1b(i,k) + stept*Dx_half(i)*DxP ) &
        /(2.0_SP + stept*(Ax_half(i)+Dx_half(i)))
      Vx(i,k)=Vx(i,k)+stept/rhox*((1.0_SP/Bx_half(i)-1)*DxP-P1/Bx_half(i))
      P1b(i,k)=2.0_SP*P1-P1b(i,k)

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

      DzP = ( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz

      ! correct v_z
      P3 =  (2.0_SP * P3a(i,k) + stept*Dz_half(k)*DzP ) &
        /(2.0_SP + stept*(Az_half(k)+Dz_half(k)))
      Vz(i,k)=Vz(i,k)+stept/rhoz*((1.0_SP/Bz_half(k)-1)*DzP-P3/Bz_half(k))
      P3a(i,k)=2.0_SP*P3-P3a(i,k)

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

      DzP = ( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz

      ! correct v_z
      P3 =  (2.0_SP * P3b(i,k) + stept*Dz_half(k)*DzP ) &
        /(2.0_SP + stept*(Az_half(k)+Dz_half(k)))
      Vz(i,k)=Vz(i,k)+stept/rhoz*((1.0_SP/Bz_half(k)-1)*DzP-P3/Bz_half(k))
      P3b(i,k)=2.0_SP*P3-P3b(i,k)

    end do
    end do
  end if ! isz2

end subroutine abs_cfs_momentum

subroutine abs_cfs_hook
  integer :: n,i,k
  real(SP) :: DxVx,DzVz
  real(SP) :: Vx1,Vx3,Vz1,Vz3
  real(SP) :: lam,miu,lam2mu

  if (isx1) then
    do k=nk1,nk2
    do i=ni1,ni1+abs_number(1,1)

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DxVx = ( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx

      ! correct Txx,Tzz
      Vx1 =  (2.0_SP * Vx1a(i,k) + stept*Dx_regu(i)*DxVx ) &
        /(2.0_SP + stept*(Ax_regu(i)+Dx_regu(i)))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Bx_regu(i)-1)*DxVx-Vx1/Bx_regu(i))
      Vx1a(i,k)=2.0_SP*Vx1-Vx1a(i,k)

      if (FreeOnTii .eqv. .true.) then
        if (freenode .and. k==nk2) then
          P(i,k)=P(i,k)-stept*lam**2/lam2mu*((1.0/Bx_regu(i)-1)*DxVx-Vx1/Bx_regu(i))
        end if
      end if

    end do
    end do
  end if ! isx1

  if (isx2) then
    do k=nk1,nk2
    do i=ni2-abs_number(1,2),ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DxVx = ( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx

      ! correct Txx,Tzz
      Vx1 =  (2.0_SP * Vx1b(i,k) + stept*Dx_regu(i)*DxVx ) &
        /(2.0_SP + stept*(Ax_regu(i)+Dx_regu(i)))
      P(i,k)=P(i,k)+stept*lam2mu*((1.0/Bx_regu(i)-1)*DxVx-Vx1/Bx_regu(i))
      Vx1b(i,k)=2.0_SP*Vx1-Vx1b(i,k)

      if (FreeOnTii .eqv. .true.) then
        if (freenode .and. k==nk2) then
          P(i,k)=P(i,k)-stept*lam**2/lam2mu*((1.0/Bx_regu(i)-1)*DxVx-Vx1/Bx_regu(i))
        end if
      end if

    end do
    end do
  end if ! isx2

  if (isz1) then
    do k=nk1,nk1+abs_number(2,1)
    do i=ni1,ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;


      DzVz = ( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz

      ! correct Txx,Tzz
      Vz3 =  (2.0_SP * Vz3a(i,k) + stept*Dz_regu(k)*DzVz ) &
        /(2.0_SP + stept*(Az_regu(k)+Dz_regu((k))))
      P(i,k)=P(i,k)+stept*lam2mu   *((1.0/Bz_regu(k)-1)*DzVz-Vz3/Bz_regu(k))
      Vz3a(i,k)=2.0_SP*Vz3-Vz3a(i,k)

    end do
    end do
  end if ! isz1

  if (isz2) then
    do k=nk2-abs_number(2,2),nk2
    do i=ni1,ni2

      lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

      DzVz = ( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz

      ! correct Txx,Tzz
      Vz3 =  (2.0_SP * Vz3b(i,k) + stept*Dz_regu(k)*DzVz ) &
        /(2.0_SP + stept*(Az_regu(k)+Dz_regu((k))))
      P(i,k)=P(i,k)+stept*lam2mu   *((1.0/Bz_regu(k)-1)*DzVz-Vz3/Bz_regu(k))
      Vz3b(i,k)=2.0_SP*Vz3-Vz3b(i,k)

    end do
    end do
  end if ! isz2

end subroutine abs_cfs_hook

!--------------------------------------------------------------------}

function cal_pml_dmax(L,Vp,Rpp,bmax) result(dmax)
  real(SP) :: L,Vp,Rpp,bmax
  real(SP) :: dmax
  !dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)
  dmax=-Vp/2.0_SP/L*log(Rpp)*(CONSPD+1.0_SP)*3
  !dmax=dmax*(CONSPD+CONSPB+1)/( (CONSPD+1.0_SP)*bmax+CONSPB )
end function cal_pml_dmax

function cal_pml_d(x,L,dmax) result(d)
  real(SP) :: x,L,dmax
  real(SP) :: d
  !real(SP),parameter :: p=2.0_SP

  !dmax=-(p+1.0_SP)*Vp/2.0_SP/L*log(Rpp)

  if (x<0) then
    d=0.0
  else
    d=dmax*(x/L)**CONSPD
  end if
end function cal_pml_d

function cal_pml_amax(fc) result(amax)
  real(SP) :: fc
  real(SP) :: amax
  amax=PI*fc
end function cal_pml_amax

function cal_pml_a(x,L,amax) result(a)
  real(SP) :: x,L,amax
  real(SP) :: a
  !amax=PI*fc
  if (x<0) then
    a=0.0
  else
    !a=amax*(L-x)/L
    a=amax*(1.0_SP-(x/L)**CONSPA)
  end if
end function cal_pml_a

function cal_pml_b(x,L,bmax) result(b)
  real(SP) :: x,L,bmax
  real(SP) :: b
  !real(SP),parameter :: p=2.0_SP
  if (x<0) then
    b=1.0
  else
    b=1.0+(bmax-1.0)*(x/L)**CONSPB
  end if
end function cal_pml_b
function cal_pml_R(n) result(r)
  integer :: n
  real(SP) :: r
  r = real(10**(-( log10(real(N,DP)) - 1.0_DP)/log10(2.0_DP) -3.0_DP),SP)
end function cal_pml_R

!---------------------- io related -----------------------------
subroutine abs_cfs_rest_export(pnm_rest)
  character (len=*) :: pnm_rest

end subroutine abs_cfs_rest_export

subroutine abs_cfs_rest_import(pnm_rest)
  character (len=*) :: pnm_rest


end subroutine abs_cfs_rest_import

end module abs_cfs_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
