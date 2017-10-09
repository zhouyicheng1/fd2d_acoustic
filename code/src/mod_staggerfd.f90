

!******************************************************************************!
!*                                                                            *!
!******************************************************************************!
!----------------------------------------------------------------------------
! This module contains the variables and subroutines used in the 
! staggered fd operator                                 
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
module solver_mod
!----------------------------------------------------------------------------
































































!#define VPointSymm

use constants_mod, only : SEIS_GEO
use para_mod
use nompi_mod
use media_mod
use grid_mod

implicit none
private
public ::                   &
  P,Vx,Vz,     &
  update_momentum,          &
  update_hook,              &
  solver_init,              &
  solver_destroy,           &
  solver_check,             &
  atten_T, atten_V

real,parameter,private :: MCC0=1.125, MCC1=0.0416666667

real(SP),dimension(:,:),allocatable :: P,Vx,Vz
real(SP),dimension(:),allocatable,public :: PSrc,VxSrc,VzSrc
integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: indx

!----------------------------------------------------------------------------
contains
!----------------------------------------------------------------------------

subroutine solver_init
  integer :: ierr
  allocate( P(nx1:nx2,nz1:nz2),stat=ierr);  P=0.0
  allocate( Vx (nx1:nx2,nz1:nz2),stat=ierr);  Vx =0.0
  allocate( Vz (nx1:nx2,nz1:nz2),stat=ierr);  Vz =0.0
  allocate(PSrc(nx1:nx2),stat=ierr); PSrc=0.0
  allocate(VxSrc(nx1:nx2),stat=ierr); VxSrc=0.0
  allocate(VzSrc(nx1:nx2),stat=ierr); VzSrc=0.0
  if (ierr>0) then
    print *, "can't allocate variable in solver_init"
    stop 1
  end if
  ! main
  indx(:,SEIS_GEO*2+1)=(/ ni1+2,ni2-2, &
    nk1+2,nk2-2 /)
  indx(:,SEIS_GEO*2  )=(/ ni1,ni2,nk2-2+1,nk2 /) ! z2
  indx(:,SEIS_GEO*2-1)=(/ ni1,ni2,nk1,nk1+2-1 /) ! z1
  indx(:,1)=(/ ni1,ni1+2-1,nk1+2,nk2-2 /) ! x1
  indx(:,2)=(/ ni2-2+1,ni2,nk1+2,nk2-2 /) ! x2
end subroutine solver_init
subroutine solver_destroy
  deallocate( P, Vx, Vz)
end subroutine solver_destroy

subroutine solver_check(ntime)
  integer ntime
  real(SP) :: V1,V3,P1,W
  integer ierr
  if (mod(ntime,1)==0) then
    V1=maxval(abs(Vx))
    V3=maxval(abs(Vz))
    P1=maxval(abs(P))
    W=max(V1,V3,P1)
    if (W>=huge(1.0)) then
      print *, "Overflow error: "
      write(*,"(i5,i3.2,1(i2.2),5(es12.5,3i5))") ntime,thisid(1),thisid(2), &
        V1, maxloc(abs(Vx)), V3, maxloc(abs(Vz)),                          &
        P1, maxloc(abs(P))
      stop 1
    end if
  end if
end subroutine solver_check

subroutine atten_T
  integer i,j,k
  real(SP) :: Qatt
  real(SP) :: Qxz
  !$OMP PARALLEL DEFAULT(shared)
  !$OMP DO PRIVATE(i,j,k)
  do k=nk1,nk2
  do i=ni1,ni2
    if (Qs(i,k)<QsINF) then
      Qatt=exp((-PI*QsF0*stept)/Qs(i,k))
      P(i,k)=P(i,k)*Qatt
    end if
  end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine atten_T

subroutine atten_V
  integer i,j,k
  real(SP) :: Qatt
  real(SP) :: Qx,Qz
  !$OMP PARALLEL DEFAULT(shared)
  !$OMP DO PRIVATE(i,j,k,Qx)
  do k=nk1,nk2
  do i=ni1,ni2
    if (Qs(i,k)<QsINF) then
      Qx=2.0/(Qs(i,k)+1.0/Qs(i+1,k))
      Qatt=exp((-PI*QsF0*stept)/Qx)
      Vx (i,k)=Vx (i,k)*Qx
    end if
  end do
  end do
  !$OMP END DO
  !$OMP DO PRIVATE(i,j,k,Qz)
  do k=nk1,nk2
  do i=ni1,ni2
    if (Qs(i,k)<QsINF) then
      Qz=2.0/(Qs(i,k)+1.0/Qs(i,k+1))
      Qatt=exp((-PI*QsF0*stept)/Qz)
      Vz (i,k)=Vz (i,k)*Qz
    end if
  end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
end subroutine atten_V

!-------------------------------------------------------
subroutine update_momentum
  integer :: n,ierr
  integer :: k0,i

  if ( freenode ) then
    if(CondFreeSIMG .eqv. .true.) then
      call planar_simg
    end if
    if(CondFreeAFDA .eqv. .true.) then
      call planar_safda
    end if
  end if

  n=SEIS_GEO*2+1
  call cal_momentum( indx(1,n),indx(2,n), &
    indx(3,n),indx(4,n) )

  do n=1,SEIS_GEO*2
  call cal_momentum( indx(1,n),indx(2,n), &
    indx(3,n),indx(4,n) )
  end do

end subroutine update_momentum

subroutine planar_simg
  integer i,j,k0

  if(FreeOnTii .eqv. .true.) then
    ! T=0
    do i=ni1,ni2
      P(i,nk2)=PSrc(i)
    end do
    ! T anti-symmetric
    !-- do k0=1,2-1
    do k0=1,2
    do i=ni1,ni2
      P(i,nk2+k0)=2.0*PSrc(i)-P(i,nk2-k0)
    end do
    end do

    ! Tzz anti-symmetric
    do k0=1,2
    do i=ni1,ni2
      P(i,nk2+k0)=-P(i,nk2-k0+1)
    end do
    end do

  end if

end subroutine planar_simg

subroutine planar_safda
  integer i,k0

  if (FreeOnTii .eqv. .true.) then

    ! Tzz=0
    do i=ni1,ni2
      P(i,nk2)=PSrc(i)
    end do


  end if
end subroutine planar_safda

subroutine cal_momentum(I1,I2,K1,K2)
  integer,intent(in) :: I1,I2,K1,K2
  integer :: i,k
  real(SP) :: rhox,rhoz


  ! ========================  Vx ========================
  ! === interior domain upto nk2-2 ===
  !$OMP PARALLEL DEFAULT(shared)
  !$OMP DO PRIVATE(i,j,k,rhox)
  do k=K1,min(K2,nk2-2)
  do i=I1,I2
    if(MediaPrecise .eqv. .true.) then
      rhox=Px(i,k)
    else
      rhox=0.5*(rho(i,k)+rho(i+1,k))
    end if
    Vx(i,k)=Vx(i,k) - stept*( &
      ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx         &
      )/rhox
  end do
  end do
  !$OMP END DO

  ! === above nk2-2  ===
  !$OMP PARALLEL DEFAULT(shared)
  !$OMP DO PRIVATE(i,j,k,rhox)
  do k=nk2-2+1,K2
  do i=I1,I2
    if(MediaPrecise .eqv. .true.) then
      rhox=Px(i,k)
    else
      rhox=0.5*(rho(i,k)+rho(i+1,k))
    end if

    if(CondFreeSIMG .eqv. .true.) then
      Vx(i,k)=Vx(i,k) - stept*( &
        ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx         &
        )/rhox
    end if

    if(CondFreeAFDA .eqv. .true.) then

      if (FreeOnTii .eqv. .true.) then

        if (freenode .and. k==nk2) then
          Vx(i,k)=Vx(i,k) - stept*( &
            ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx         &
           )/rhox
        elseif (freenode .and. k==nk2-1) then
          Vx(i,k)=Vx(i,k) - stept*( &
            ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx         &
            )/rhox
        else
          Vx(i,k)=Vx(i,k) - stept*( &
            ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx         &
            )/rhox
      end if

      else

        if (freenode .and. k==nk2) then
          Vx(i,k)=Vx(i,k) - stept*( &
            ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx       &
           )/rhox
        else
          Vx(i,k)=Vx(i,k) - stept*( &
            ( MCC0*(P(i+1,k)-P(i,k)) -MCC1*(P(i+2,k)-P(i-1,k)) )/stepx       &
            )/rhox
        end if

      end if

    end if

  end do
  end do
  !$OMP END DO

  ! ========================  Vz ========================
  !$OMP DO PRIVATE(i,j,k,rhoz)
  do k=K1,min(K2,nk2-2)
  do i=I1,I2
    if (MediaPrecise .eqv. .true.) then
      rhoz=Pz(i,k)
    else
      rhoz=0.5*(rho(i,k)+rho(i,k+1))
    end if
    Vz(i,k)=Vz(i,k) - stept*( &
      +( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz         &
      )/rhoz
  end do
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(i,j,k,rhoz)
  do k=nk2-2+1,K2
  do i=I1,I2
    if (MediaPrecise .eqv. .true.) then
      rhoz=Pz(i,k)
    else
      rhoz=0.5*(rho(i,k)+rho(i,k+1))
    end if

    if (CondFreeSIMG .eqv. .true.) then
      Vz(i,k)=Vz(i,k) - stept*( &
        +( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz         &
        )/rhoz
    end if

    if (CondFreeAFDA .eqv. .true.) then

      if(FreeOnTii .eqv. .true.) then

        if (freenode .and. k==nk2) then
          Vz(i,k)=0.0_SP
        elseif (freenode .and. k==nk2-1) then
          Vz(i,k)=Vz(i,k) - stept*( &
            +( 11./12.*P(i,k+1) -17./24.*P(i,k)-3./8.*P(i,k-1)+5./24.*P(i,k-2)-1./24.*P(i,k-3) )/stepz    &
            )/rhoz
        else
          Vz(i,k)=Vz(i,k) - stept*( &
            +( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz         &
            )/rhoz
        end if

      else

        if (freenode .and. k==nk2) then
        Vz(i,k)=Vz(i,k) - stept*( &
          +( -35./8.*P(i,k) +35./24.*P(i,k-1)-21./40.*P(i,k-2)+5./56.*P(i,k-3) )/stepz   &
          )/rhoz
        elseif (freenode .and. k==nk2-1) then
          Vz(i,k)=Vz(i,k) - stept*( &
            +( 31./24.*P(i,k+1) -29./24.*P(i,k)+3./40.*P(i,k-1)-1./168.*P(i,k-2) )/stepz   &
            )/rhoz
        else
          Vz(i,k)=Vz(i,k) - stept*( &
            +( MCC0*(P(i,k+1)-P(i,k)) -MCC1*(P(i,k+2)-P(i,k-1)) )/stepz         &
            )/rhoz
       end if
    end if
    end if

  end do
  end do
  !print*,Vz(150,135)
  !$OMP END DO

  !$OMP END PARALLEL
end subroutine cal_momentum
!-------------------------------------------------------

!-------------------------------------------------------
subroutine update_hook
  integer :: n,ierr
  integer :: k0,i

  if ( freenode ) then
    if (CondFreeV2TH .eqv. .true.) then
      if(FreeOnTii .eqv. .true.) then
        call apply_v2th_OnTii
      else
        call apply_v2th_OnTij
      end if
    end if
    if (CondFreeVZERO .eqv. .true.) then
      call apply_vzero
    end if
    if (CondFreeVEXT .eqv. .true.) then
      call swmpi_except("free surface condition vext hasn't been implemented")
      call apply_vext
    end if
  end if

  n=SEIS_GEO*2+1
  call cal_hook( indx(1,n),indx(2,n),  &
    indx(3,n),indx(4,n) )

  do n=1,SEIS_GEO*2
    call cal_hook( indx(1,n),indx(2,n),  &
      indx(3,n),indx(4,n) )
  end do

end subroutine update_hook

subroutine cal_hook(I1,I2,K1,K2)
  integer,intent(in) :: I1,I2,K1,K2
  real(SP) :: lam,miu,lam2mu
  integer i,k
  real(SP) :: DxVx,DzVz
  real(SP) :: DzVx    ! for AFDA
  integer :: k0

  !$OMP PARALLEL DEFAULT(shared)

  ! ========================  Tii ========================
  !$OMP DO PRIVATE(i,j,k,lam,miu,lam2mu)
  do k=K1,min(K2,nk2-2)
  do i=I1,I2
    lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;
    DxVx=( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx
    DzVz=( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz
    P(i,k)=P(i,k) - stept*      &
     lam2mu* (DxVx +  DzVz )
  end do
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(i,j,k,lam,miu,lam2mu)
  do k=nk2-2+1,K2
  do i=I1,I2
    lam=lambda(i,k);miu=mu(i,k);lam2mu=lam+2.0*miu;

    DxVx=( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx

    if(FreeOnTii .eqv. .true.) then

      DzVz=( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz

      if (freenode .and. k==nk2) then
        DzVz=-VzSrc(i)/lam2mu-lam/lam2mu* DxVx
      end if

      if (CondFreeAFDA .eqv. .true.) then

        if (freenode .and. k==nk2-1) then
          ! derivative at nk2
          DxVx=( MCC0*(Vx(i,k+1)-Vx(i-1,k+1)) -MCC1*(Vx(i+1,k+1)-Vx(i-2,k+1)) )/stepx
          DzVz=VzSrc(i)/lam2mu-lam/lam2mu* DxVx
          ! using compact scheme with the derivative at nk2 
          DzVz=( -stepz/22.*DzVz +577./528.*Vz(i,k) -201./176.*Vz(i,k-1)+9./176.*Vz(i,k-2)-1./528.*Vz(i,k-3) )/stepz
          DxVx=( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx
        end if

      end if

    else

      if (freenode .and. k==nk2) then

        if (CondFreeAFDA .eqv. .true.) then
          DzVz=( 11./12.*Vz(i,k) -17./24.*Vz(i,k-1)-3./8.*Vz(i,k-2)+5./24.*Vz(i,k-3)-1./24.*Vz(i,k-4) )/stepz
        else
          DzVz=( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz
        end if
      else
        DzVz=( MCC0*(Vz(i,k)-Vz(i,k-1)) -MCC1*(Vz(i,k+1)-Vz(i,k-2)) )/stepz
      end if
    end if

    P(i,k)=P(i,k) - stept*      &
      lam2mu* (DxVx + DzVz )
  end do
  end do
  !$OMP END DO

  ! ========================  Txz ========================
  !$OMP DO PRIVATE(i,j,k,miu)
  do k=K1,min(K2,nk2-2)
  do i=I1,I2
    if (MediaPrecise .eqv. .true.) then
      miu=Uxz(i,k)
    else

      if (mu(i,k)<=SEIS_ZERO .or. mu(i,k+1)<=SEIS_ZERO .or. mu(i+1,k)<=SEIS_ZERO .or. mu(i+1,k+1)<=SEIS_ZERO) then
        miu=0
      else
        miu=4.0/(1.0/mu(i,k)+1.0/mu(i,k+1)+1.0/mu(i+1,k)+1.0/mu(i+1,k+1))
      end if

    end if
    !-- #endif
!    Txz(i,k)=Txz(i,k) + stept*miu*( &
!      ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx             &
!      +( MCC0*(Vx(i,k+1)-Vx(i,k)) -MCC1*(Vx(i,k+2)-Vx(i,k-1)) )/stepz )
  end do
  end do
  !$OMP END DO

  !$OMP DO PRIVATE(i,j,k,miu)
  do k=nk2-2+1,K2
  do i=I1,I2
    if(MediaPrecise .eqv. .true.) then
      miu=Uxz(i,k)
    else
      if (mu(i,k)<=SEIS_ZERO .or. mu(i,k+1)<=SEIS_ZERO .or. mu(i+1,k)<=SEIS_ZERO .or. mu(i+1,k+1)<=SEIS_ZERO) then
        miu=0
      else
        miu=4.0/(1.0/mu(i,k)+1.0/mu(i,k+1)+1.0/mu(i+1,k)+1.0/mu(i+1,k+1))
      end if
      !--miu=4.0/(1.0/mu(i,k)+1.0/mu(i,k+1)+1.0/mu(i+1,k)+1.0/mu(i+1,k+1))
    end if

    if (CondFreeAFDA .eqv. .true.) then

      if(FreeOnTii .eqv. .true.) then

        if (freenode .and. k==nk2) then
 !         Txz(i,k)=0.0_SP
        elseif (freenode .and. k==nk2-1) then
 !         Txz(i,k)=Txz(i,k) + stept*miu*( &
 !           ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx                &
 !           +( 11./12.*Vx(i,k+1) -17./24.*Vx(i,k)-3./8.*Vx(i,k-1)+5./24.*Vx(i,k-2)-1./24.*Vx(i,k-3) )/stepz )
        else
 !         Txz(i,k)=Txz(i,k) + stept*miu*( &
 !           ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx                &
 !           +( MCC0*(Vx(i,k+1)-Vx(i,k)) -MCC1*(Vx(i,k+2)-Vx(i,k-1)) )/stepz )
        end if

      else

        if (freenode .and. k==nk2) then
 !         Txz(i,k)=0.0_SP
        elseif (freenode .and. k==nk2-1) then
          DzVx=-( MCC0*(Vz(i+1,k+1)-Vz(i,k+1)) -MCC1*(Vz(i+2,k+1)-Vz(i-1,k+1)) )/stepx
 !         Txz(i,k)=Txz(i,k) + stept*miu*( &
 !           ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx                &
 !           +( -stepz/22.*DzVx +577./528.*Vx(i,k+1) -201./176.*Vx(i,k)+9./176.*Vx(i,k-1)-1./528.*Vx(i,k-2) )/stepz )
        else
 !         Txz(i,k)=Txz(i,k) + stept*miu*( &
 !           ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx                &
 !           +( MCC0*(Vx(i,k+1)-Vx(i,k)) -MCC1*(Vx(i,k+2)-Vx(i,k-1)) )/stepz )
        end if
        
      end if

    else
 !     Txz(i,k)=Txz(i,k) + stept*miu*(  &
 !       ( MCC0*(Vz(i+1,k)-Vz(i,k)) -MCC1*(Vz(i+2,k)-Vz(i-1,k)) )/stepx                &
 !       +( MCC0*(Vx(i,k+1)-Vx(i,k)) -MCC1*(Vx(i,k+2)-Vx(i,k-1)) )/stepz )
    end if

  end do
  end do
  !$OMP END DO

  !$OMP END PARALLEL
end subroutine cal_hook

subroutine symmetric_cond
  integer i
  do i=ni1,ni2
    Vx (i,nk2+1:nk2+2)=Vx (i,nk2-1:nk2-2:-1)
    Vz (i,nk2+1:nk2+2)=Vz (i,nk2-1:nk2-2:-1)
  end do
end subroutine symmetric_cond

subroutine apply_vzero
  integer i
  do i=ni1,ni2
    Vx (i,nk2+1:nk2+2)=0.0
    Vz (i,nk2+1:nk2+2)=0.0
    if(FreeOnTii .eqv. .true.) then
      Vz (i,nk2:nk2+2)=0.0
    end if
  end do
end subroutine apply_vzero

subroutine apply_v2th_OnTii
  real(SP),parameter :: c0=1.125, c1=0.0416666667
  integer i,k
  real(SP) :: dz,dx,lam,miu,lam2mu
  !Vz
  do i=ni1,ni2
    k=nk2
    dx=(x(i+1)-x(i-1))/2.0
    dz=(z(k+1)-z(k-1))/2.0

    lam=lambda(i,k); miu=mu(i,k);lam2mu=lam+2.0*miu;
    Vz(i,k)=-dz*lam/lam2mu*(Vx(i,k)-Vx(i-1,k))/dx &
      +dz/lam2mu*VzSrc(i)                     &
      +Vz(i,k-1)
    Vz(i,k+1)=(c0*(Vz(i,k)-Vz(i,k-1))             &
      +dz*lam/lam2mu*( MCC0*(Vx(i,k)-Vx(i-1,k)) -MCC1*(Vx(i+1,k)-Vx(i-2,k)) )/stepx      &
      -dz/lam2mu*VzSrc(i))/c1              &
      +Vz(i,k-2)
  end do
  !Vx
  do i=ni1,ni2
    k=nk2
    miu=mu(i,k)
    dz=z(k+1)-z(k)
    dx=x(i+1)-x(i)
    Vx(i,k+1)=( -(Vx(i,k)-Vx(i,k-1))/dz-(Vz(i+1,k-1)-Vz(i,k-1))/dx &
      +2.0*VxSrc(i)/miu                                         &
      -(Vz(i+1,k)-Vz(i,k))/dx)*dz                            &
      +Vx(i,k)
  end do
end subroutine apply_v2th_OnTii
subroutine apply_v2th_OnTij
  integer i,k
  real(SP) :: dz,dx,lam,miu
  !Vx
  do i=ni1,ni2
    k=nk2
    dz=z(k+1)-z(k)
    dx=x(i+1)-x(i)
    Vx(i,k+1)=Vx(i,k)-dz/dx*(Vz(i+1,k)-Vz(i,k))
  end do
  !Vz
  do i=ni1,ni2
    dx=(x(i+1)-x(i-1))/2.0

    k=nk2
    dz=(z(k+1)-z(k-1))/2.0
    lam=lambda(i,k); miu=mu(i,k)
    Vz(i,k+1)=-lam/(lam+2.0*miu)*(Vx(i,k)-Vx(i-1,k))/dx &
      -(Vz(i,k)-Vz(i,k-1))/dz

    k=nk2+1
    dz=(z(k+1)-z(k-1))/2.0
    Vz(i,k)=Vz(i,k)*dz+Vz(i,k-1) -lam/(lam+2.0*miu) &
      *( dz/dx*(Vx(i,k)-Vx(i-1,k)) )
  end do
end subroutine apply_v2th_OnTij

subroutine apply_vext
  return
end subroutine

!--------------------------------------------------------------------}

end module solver_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
