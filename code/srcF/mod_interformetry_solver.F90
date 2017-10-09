module interformetry_solver_mod
!
!
! This module is used for interformetry
!
! Wei Zhang     Email: zhangwei.zw@gmail.con
! Yicheng Zhou  Email: zhouyc07@mail.ustc.edu.cn
! Copyright(C) 2006 Wei ZHANG          
!
!*********************************************************
!
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $      
! $Date: 2009-01-14 16:30:51 -0500 (Wed, 14 Jan 2009) $      
! $Revision: 510 $       
! $LastChangeBy: yicheng zhou $
!
!*********************************************************
!
use constants_mod
use string_mod
use math_mod
use para_mod
use nompi_mod
use grid_mod
use media_mod
use nfseis_mod
use interformetry_mod
use solver_mod,  only: Vpx,Vpz,Vsx,Vsz
implicit none
private

public ::             &
  SVpx,SVpz,SVsx,SVsz,&
  GVpx,GVpz,GVsx,GVsz,&
  IVpx,IVpz,IVsx,IVsz,&
  IGVpx,IGVpz,IGVsx,IGVsz,&
  Imgcon,IImgcon,&
  solver_interformetry_get, &
  solver_interformetry_init,&
  update_interformetry,&
  update_img
real(SP),dimension(:,:),allocatable :: SVpx, SVpz,SVsx,SVsz,GVpx,GVpz,GVsx,GVsz,&
                                       IVpx,IVpz,IVsx,IVsz,IGVpx,IGVpz,IGVsx,IGVsz,&
                                       Imgcon,IImgcon
integer,dimension(SEIS_GEO*2,SEIS_GEO*2+1),public :: iindx

!---------------------------------------------------------
contains
!---------------------------------------------------------
subroutine solver_interformetry_init
  integer :: ierr
  allocate( SVpx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   SVpx=0.0
  allocate( SVpz(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   SVpz=0.0
  allocate( SVsx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   SVsx=0.0
  allocate( SVsz(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   SVsz=0.0
  allocate( IVpx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   IVpx=0.0
  allocate( IVsx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);   IVsx=0.0
  allocate( IGVpx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  IGVpx=0.0
  allocate( IGVpz(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  IGVpz=0.0
  allocate( IGVsx(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  IGVsx=0.0
  allocate( IGVsz(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  IGVsz=0.0
  allocate( Imgcon(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  Imgcon=0.0
  allocate( IImgcon(inx1+iffo1:inx2+iffo2,inz1+iffo1:inz2+iffo2),stat=ierr);  IImgcon=0.0
  if (ierr>0) then
    print *, "can't allocate variable in solver_init"
    stop 1
  end if
  ! main
!  iindx(:,SEIS_GEO*2+1)=(/ inx1+iffo2,inx2+iffo1, &
!    inz1+iffo2,inz2+iffo1 /)
!  iindx(:,SEIS_GEO*2  )=(/ inx1,inx2,inz2+iffo1+1,inz2 /) ! z2
!  iindx(:,SEIS_GEO*2-1)=(/ inx1,inx2,inz1,inz1+iffo2-1 /) ! z1
!  iindx(:,1)=(/ inx1,inx1+iffo2-1,inz1+iffo2,inz2+iffo1 /) ! x1
!  iindx(:,2)=(/ inx2+iffo1+1,inx2,inz1+iffo2,inz2+iffo1 /) ! x2
  ! main
  iindx(:,SEIS_GEO*2+1)=(/ inx1,inx2, inz1,inz2 /)
  iindx(:,SEIS_GEO*2  )=(/ inx1+iffo1,inx2+iffo2,inz2+1,inz2+iffo2 /) ! z2
  iindx(:,SEIS_GEO*2-1)=(/ inx1+iffo1,inx2+iffo2,inz1+iffo1,inz1 /) ! z1
  iindx(:,1)=(/ inx1+iffo1,inx1-1,inz1,inz2+iffo1 /) ! x1
  iindx(:,2)=(/ inx2+1,inx2+iffo2,inz1,inz2+iffo1 /) ! x2
end subroutine solver_interformetry_init

subroutine solver_interformetry_get
    integer fid
    fid=1031
    call string_conf(fid,1,'INT_ROOT',2,pnm_int)       ! pnm_int--INT_ROOT = ./input
    fnm_int=trim(pnm_int)//'/'//'vel_001_n00001.nc'
    call nfseis_varget(fnm_int,'Vpx',Vpx(:,:),(/nx1,nz1/),(/nx2,nz2/),(/1,1/))
    call nfseis_varget(fnm_int,'Vpz',Vpz(:,:),(/nx1,nz1/),(/nx2,nz2/),(/1,1/))
    call nfseis_varget(fnm_int,'Vsx',Vsx(:,:),(/nx1,nz1/),(/nx2,nz2/),(/1,1/))
    call nfseis_varget(fnm_int,'Vsz',Vsz(:,:),(/nx1,nz1/),(/nx2,nz2/),(/1,1/))
end subroutine solver_interformetry_get

subroutine update_interformetry
integer :: n,ierr
  n=SEIS_GEO*2+1
  call cal_interformetry(iffo1,iffo2, &
    iindx(1,n),iindx(2,n),  &
    iindx(3,n),iindx(4,n) )

  do n=1,SEIS_GEO*2
    call cal_interformetry(iffo1,iffo2,&
      iindx(1,n),iindx(2,n),  &
      iindx(3,n),iindx(4,n) )
  end do
end subroutine update_interformetry

subroutine update_img
  integer :: n,ierr
  n=SEIS_GEO*2+1
  do n=1,SEIS_GEO*2
     call cal_img( iindx(1,n),iindx(2,n),  &
                   iindx(3,n),iindx(4,n) )
  end do
end subroutine update_img

subroutine cal_interformetry(iffo1,iffo2,I1,I2,K1,K2)
  integer,intent(in) :: iffo1,iffo2,I1,I2,K1,K2
  integer :: i,k,ifo1,ifo2
  do k=K1,K2
  do i=I1,I2
     SVpx(i,k)=0.0
     SVsx(i,k)=0.0
     SVpz(i,k)=0.0
     SVsz(i,k)=0.0
     do ifo1=iffo1,iffo2
     do ifo2=iffo1,iffo2
        SVpx(i,k) = SVpx(i,k) + Vpx(i+ifo1,k+ifo2)
        SVsx(i,k) = SVsx(i,k) + Vsx(i+ifo1,k+ifo2)
        SVpz(i,k) = SVpz(i,k) + Vpz(i+ifo1,k+ifo2)
        SVsz(i,k) = SVsz(i,k) + Vsz(i+ifo1,k+ifo2)
     end do
     end do
  end do
  end do
  do k=K1,K2
  do i=I1,I2
        SVpx(i,k) = SVpx(i,k)/((iffo2-iffo1+1)**2)
        SVsx(i,k) = SVsx(i,k)/((iffo2-iffo1+1)**2)
        SVpz(i,k) = SVpz(i,k)/((iffo2-iffo1+1)**2)
        SVsz(i,k) = SVsz(i,k)/((iffo2-iffo1+1)**2)
  end do
  end do
  do k=K1,K2
  do i=I1,I2
        GVpx(i,k) = Vpx(i,k) - SVpx(i,k)
        GVsx(i,k) = Vsx(i,k) - SVsx(i,k)
        GVpz(i,k) = Vpz(i,k) - SVpz(i,k)
        GVsz(i,k) = Vsz(i,k) - SVsz(i,k)
  end do
  end do
  do k=K1,K2
  do i=I1,I2
     IVpx(i,k)=0.0
     IVpz(i,k)=0.0
     IVsx(i,k)=0.0
     IVsz(i,k)=0.0
     do ifo1=iffo1,iffo2
     do ifo2=iffo1,iffo2
        IVpx(i,k) = IVpx(i,k) + Vpx(i+ifo1,k+ifo2)*Vpx(i-ifo1,k-ifo2)
        IVsx(i,k) = IVsx(i,k) + Vsx(i+ifo1,k+ifo2)*Vsx(i-ifo1,k-ifo2)
        IVpz(i,k) = IVpz(i,k) + Vpz(i+ifo1,k+ifo2)*Vpz(i-ifo1,k-ifo2)
        IVsz(i,k) = IVsz(i,k) + Vsz(i+ifo1,k+ifo2)*Vsz(i-ifo1,k-ifo2)
     end do
     end do
  end do
  end do
  do k=K1,K2
  do i=I1,I2
     IGVpx(i,k)=0.0
     IGVpz(i,k)=0.0
     IGVsx(i,k)=0.0
     IGVsz(i,k)=0.0
     do ifo1=iffo1,iffo2
     do ifo2=iffo1,iffo2
        IGVpx(i,k) = IGVpx(i,k) + GVpx(i+ifo1,k+ifo2)*GVpx(i-ifo1,k-ifo2)
        IGVsx(i,k) = IGVsx(i,k) + GVsx(i+ifo1,k+ifo2)*GVsx(i-ifo1,k-ifo2)
        IGVpz(i,k) = IGVpz(i,k) + GVpz(i+ifo1,k+ifo2)*GVpz(i-ifo1,k-ifo2)
        IGVsz(i,k) = IGVsz(i,k) + GVsz(i+ifo1,k+ifo2)*GVsz(i-ifo1,k-ifo2)
     end do
     end do
  end do
  end do
end subroutine cal_interformetry

subroutine cal_img(I1,I2,K1,K2)
  integer,intent(in) :: I1,I2,K1,K2
  integer :: i,k
  do i=I1,I2
  do k=K1,K2
     Imgcon(i,k)=(((Vpx(i-1,k)+Vpx(i,k))/2)**2+((Vpz(i,k-1)+Vpz(i,k))/2)**2)*&
            (((Vsx(i-1,k)+Vsx(i,k))/2)**2+((Vsz(i,k-1)+Vsz(i,k))/2)**2)
     IImgcon(i,k)=(((IVpx(i-1,k)+IVpx(i,k))/2)**2+((IVpz(i,k-1)+IVpz(i,k))/2)**2)*&
             (((IVsx(i-1,k)+IVsx(i,k))/2)**2+((IVsz(i,k-1)+IVsz(i,k))/2)**2)
  end do
  end do
end subroutine cal_img

end module interformetry_solver_mod
