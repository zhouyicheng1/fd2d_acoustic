module io_interformetry_mod

!----------------------------------------------------------------------------
! This module is used for input and output data
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
use nompi_mod
use interformetry_solver_mod
implicit none
private
public :: io_interformetry_snap_read,io_wave_interformetry_export,&
          io_interformetry_wave_close
#ifdef IOAuxilaryVar
real(SP),dimension(:,:),allocatable :: iiovar
#endif
integer,allocatable,public :: &
  snap_int_tinv(:),               &
  snap_int_tcnt(:),               &
  snap_int_gsubs(:,:),            &
  snap_int_gsubt(:,:),            &
  snap_int_gsube(:,:),            &
  snap_int_gsubc(:,:),            &
  snap_int_subs(:,:),             &
  snap_int_subt(:,:),             &
  snap_int_sube(:,:),             &
  snap_int_subc(:,:),             &
  sgt_int_ncid(:),                &
  sgt_int_vid(:,:),               &
  sgt_int_tid(:),                 &
  vel_int_ncid(:),                &
  vel_int_vid(:,:),               &
  vel_int_tid(:)
logical,allocatable,public :: &
  snap_int_ishere(:),             &
  snap_int_onfree(:),             &
  sgt_int_oflag(:),               &
  vel_int_oflag(:),               &
  sgt_int_out(:),                 &
  vel_int_out(:)

subroutine io_snap_read(fnm_int_input)
  character (len=*),intent(in) :: fnm_int_input
  character (len=2) :: varnm

  open(fid,file=trim(fnm_input),status="old")
  !-- output --
  call string_conf(fid,1,'number_of_int_snap',2,num_int_snap)
  call alloc_int_snap(num_int_snap)
  do n=1,num_int_snap
    do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('snap_',n)),m+1,snap_subs(m,n))
      call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+SEIS_GEO,snap_subc(m,n))
      call string_conf(fid,1,trim(io_enum('snap_',n)),m+1+2*SEIS_GEO,snap_subt(m,n))
    end do
    call string_conf(fid,1,trim(io_enum('snap_',n)),1+1+3*SEIS_GEO,snap_tinv(n))
    call string_conf(fid,1,trim(io_enum('snap_',n)),2+1+3*SEIS_GEO,snap_tcnt(n))
    call string_conf(fid,1,trim(io_enum('snap_',n)),3+1+3*SEIS_GEO,varnm)
    sgt_out(n)=.false.; if (index(varnm,"T")/=0) sgt_out(n)=.true.
    vel_out(n)=.false.; if (index(varnm,"V")/=0) vel_out(n)=.true.
    call corr_indx(snap_subs(:,n),snap_subc(:,n),snap_subt(:,n),snap_sube(:,n))
  end do
  !snap_sube=snap_subs+snap_subt*(snap_subc-1)
  close(fid)
end subroutine io_snap_read

subroutine alloc_int_snap(nsnap)
  integer,intent(in) :: nsnap
  allocate(snap_int_gsubs(SEIS_GEO,nsnap))
  allocate(snap_int_gsube(SEIS_GEO,nsnap))
  allocate(snap_int_gsubt(SEIS_GEO,nsnap))
  allocate(snap_int_gsubc(SEIS_GEO,nsnap))
  allocate(snap_int_subs(SEIS_GEO,nsnap))
  allocate(snap_int_sube(SEIS_GEO,nsnap))
  allocate(snap_int_subt(SEIS_GEO,nsnap))
  allocate(snap_int_subc(SEIS_GEO,nsnap))
  allocate(snap_int_tinv(nsnap))
  allocate(snap_int_tcnt(nsnap))
  allocate(snap_int_ishere(nsnap)); snap_int_ishere=.false.
  allocate(snap_int_onfree(nsnap)); snap_int_onfree=.false.
  allocate(vel_int_oflag(nsnap)); vel_int_oflag=.false.
  allocate(vel_int_ncid(nsnap))
  allocate(vel_int_vid(SEIS_GEO,nsnap))
  allocate(vel_int_tid(nsnap))
  allocate(vel_int_out(nsnap)); vel_int_out=.false.
  allocate(sgt_int_oflag(nsnap)); sgt_int_oflag=.false.
  allocate(sgt_int_ncid(nsnap))
  allocate(sgt_int_vid(2*SEIS_GEO,nsnap))
  allocate(sgt_int_tid(nsnap))
  allocate(sgt_int_out(nsnap)); sgt_int_out=.false.
end subroutine alloc_int_snap
subroutine io_interformetry_snap_locate(n_i,n_k)
  integer,intent(in) :: n_i,n_K
  integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
  integer,dimension(SEIS_GEO) :: gsubs,gsube
  integer n

  do n=1,num_snap
    subs=snap_int_subs(:,n);subc=snap_int_subc(:,n);subt=snap_int_subt(:,n)
    sube=snap_int_sube(:,n)
    ! convert into this thread
    if (      subs(1)<=ngi2 .and. sube(1)>=ngi1 &
      .and. subs(2)<=ngk2 .and. sube(2)>=ngk1 ) then
      snap_int_ishere(n)=.true.
      if (freenode .and. subs(2)==ngk2) snap_onfree(n)=.true.
      call corr_subse(subs,subc,subt,sube)
      gsubs(1)=out_i(subs(1)); gsubs(2)=out_k(subs(2))
      gsube(1)=out_i(sube(1)); gsube(2)=out_k(sube(2))
      snap_int_gsubs(:,n)=gsubs
      snap_int_gsube(:,n)=gsube
      snap_int_gsubc(:,n)=subc
      snap_int_subs(1,n)=swmpi_locli(subs(1),n_i)
      snap_int_subs(2,n)=swmpi_loclk(subs(2),n_k)
      snap_int_sube(1,n)=swmpi_locli(sube(1),n_i)
      snap_int_sube(2,n)=swmpi_loclk(sube(2),n_k)
      snap_int_subc(:,n)=subc
    end if
  end do
  call alloc_V0
#ifdef IOAuxilaryVar
  allocate(iiovar(inx1:inx2,inz1:inz2)); iiovar=0.0
#endif
end subroutine io_snap_locate
subroutine wave_interformetr_export(SVpx,SVpz,SVsx,SVsz,GVpx,GVpz,GVsx,GVsz,&
                                    IVpx,IVpz,IVsx,IVsz,IGVpx,IGVpz,IGVsx,IGVsz,&
                                    ntime)
  real(SP),dimension(:,:),intent(in) :: SVpx,SVpz,SVsx,SVsz,GVpx,GVpz,GVsx,GVsz,&
                                        IVpx,IVpz,IVsx,IVsz,IGVpx,IGVpz,IGVsx,IGVsz
  integer,intent(in) :: ntime
  integer,dimension(SEIS_GEO+1) :: startput,countput,strideput
  real(SP),intent(in) :: stept

  integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
  integer,dimension(SEIS_GEO) :: isubs,isube
  character (len=SEIS_STRLEN) :: filenm
  real(SP) :: t
  integer n,n1
  #ifdef IOAuxilaryVar
     integer i,k,ierr,si,sk
     real(SP) :: p(1)
  #endif
  t=(ntime-0.5)*stept
  do n=1,num_snap
    n1=mod(ntime/snap_int_tinv(n)-1,snap_int_tcnt(n))
    subs=snap_int_subs(:,n); subt=snap_int_subt(:,n); sube=snap_int_sube(:,n)
    subc=snap_int_subc(:,n)
    isubs(1) =out_i(subs(1));isubs(2)=out_k(subs(2))
    isube(1) =out_i(sube(1));isube(2)=out_k(sube(2))

  filenm=get_fnm_snapnode(pnm_out,'vel_int_',n,ntime,thisid(1),thisid(2))
  call nfseis_snap_def(filenm,vel_int_ncid(n),vel_int_tid(n),&
                      stept*snap_tinv(n), subc,"snap of interformetric velocity feilds")
  call nfseis_snap_attdef(vel_ncid(n), &
       snap_int_gsubs(:,n), snap_int_gsubc(:,n), snap_int_gsubt(:,n), snap_int_gsube(:,n), &
       isubs, subc, subt, isube)
  call nfseis_snap_defvar(vel_int_ncid(n),'SVpx',vel_int_vid(1,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'SVpz',vel_int_vid(2,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'SVsx',vel_int_vid(3,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'SVsz',vel_int_vid(4,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'GVpx',vel_int_vid(5,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'GVpz',vel_int_vid(6,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'GVsx',vel_int_vid(7,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'GVsz',vel_int_vid(8,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IVpx',vel_int_vid(9,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IVpz',vel_int_vid(10,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IVsx',vel_int_vid(11,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IVsz',vel_int_vid(12,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IGVpx',vel_int_vid(13,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IGVpz',vel_int_vid(14,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IGVsx',vel_int_vid(15,n))
  call nfseis_snap_defvar(vel_int_ncid(n),'IGVsz',vel_int_vid(16,n))
  call nfseis_snap_enddef(vel_int_ncid(n))
  vel_int_oflag(n)=.true.

    startput=(/1,1,n1+1/);countput=(/subc,1/);strideput=(/1,1,1/)

    call nfseis_put(vel_int_ncid(n),vel_int_tid(n),t,(/n1+1/),(/1/),(/1/) )

    if (ResultOnColocatedGrid .eqv. .true.) then
#ifdef IOAuxilaryVar
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(SVpx(i,k)+SVpx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(1,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(SVpz(i,k)+SVpz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(2,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(SVsx(i,k)+SVsx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(3,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(SVsz(i,k)+SVsz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(4,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(GVpx(i,k)+GVpx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(5,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(GVpz(i,k)+GVpz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(6,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(GVsx(i,k)+GVsx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(7,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(GVsz(i,k)+GVsz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(8,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IVpx(i,k)+IVpx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(9,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IVpz(i,k)+IVpz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(10,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IVsx(i,k)+IVsx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(11,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IVsz(i,k)+IVsz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(12,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IGVpx(i,k)+IGVpx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(13,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IGVpz(i,k)+IGVpz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(14,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IGVsx(i,k)+IGVsx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(15,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(IGVsz(i,k)+IGVsz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(16,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
#else
      sk=0
      do k=subs(2),sube(2),subt(2)
        sk=sk+1; si=0
        do i=subs(1),sube(1),subt(1)
          si=si+1
          p(1)=(SVpx(i,k)+SVpx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(1,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(SVpz(i,k)+SVpz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(2,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(SVsx(i,k)+SVsx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(3,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(SVsz(i,k)+SVsz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(4,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(GVpx(i,k)+GVpx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(5,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(GVpz(i,k)+GVpz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(6,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(GVsx(i,k)+GVsx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(7,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(GVsz(i,k)+GVsz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(8,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IVpx(i,k)+IVpx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(9,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IVpz(i,k)+IVpz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(10,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IVsx(i,k)+IVsx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(11,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IVsz(i,k)+IVsz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(12,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IGVpx(i,k)+IGVpx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(13,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IGVpz(i,k)+IGVpz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(14,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IGVsx(i,k)+IGVsx(i-1,k))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(15,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
          p(1)=(IGVsz(i,k)+IGVsz(i,k-1))/2.0
          ierr=nf90_put_var(vel_int_ncid(n),&
          vel_int_vid(16,n), p,   &
            (/si,sk,n1+1/),(/1,1,1/),(/1,1,1/) )
        end do
      end do
#endif
    else
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(1,n), &
        SVpx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(2,n), &
        SVpz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(3,n), &
        SVsx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(4,n), &
        SVsz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(5,n), &
        GVpx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(6,n), &
        GVpz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(7,n), &
        GVsx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(8,n), &
        GVsz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(9,n), &
        IVpx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(10,n), &
        IVpz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(11,n), &
        IVsx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(12,n), &
        IVsz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(13,n), &
        IGVpx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(14,n), &
        IGVpz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(15,n), &
        IGVsx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_int_ncid(n),&
        vel_int_vid(16,n), &
        IGVsz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
    end if
    if ( n1==snap_int_tcnt(n)-1 ) then
      call nfseis_close(vel_int_ncid(n))
      vel_int_oflag(n)=.false.
    end if

  end do !n

end subroutine wave_interformetry_export
subroutine io_interformetry_wave_close
  integer :: n
  do n=1,num_snap
    if (vel_int_oflag(n)) call nfseis_close(vel_int_ncid(n))
  end do
end subroutine io_interformetry_wave_close

end module io_interformetry_mod
