

module io_mod

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
implicit none

private
public :: io_init
public :: io_pt_read,io_pt_import,io_snap_read,io_snap_locate
public :: io_rest_export,io_rest_import
public :: io_seismo_init,io_seismo_put,io_seismo_close
public :: io_wave_export,io_wave_close

public :: io_stress_export,io_stress_close

public :: io_enum,io_out_pattern
public :: io_delete
public :: get_fnm_seismo,get_fnm_station,   &
  get_fnm_snapnode,get_fnm_snapnode_n
public :: corr_subse,corr_indx

character (len=SEIS_STRLEN),public ::       &
  fnm_log

real(SP),dimension(:,:),allocatable :: iovar


!---------------------------------------------
integer,public :: num_recv,num_line,num_snap,num_pt

integer,public :: pt_tinv
integer,allocatable,public :: &
  pt_indx(:,:),               &
  pt_id  (:,:),               &
  num_line_pt(:)
real(SP),allocatable,public :: pt_xz(:,:)
real(SP),public :: topo_hyper_height

integer,allocatable,public :: &
  snap_tinv(:),               &
  snap_tcnt(:),               &
  snap_gsubs(:,:),            &
  snap_gsubt(:,:),            &
  snap_gsube(:,:),            &
  snap_gsubc(:,:),            &
  snap_subs(:,:),             &
  snap_subt(:,:),             &
  snap_sube(:,:),             &
  snap_subc(:,:),             &
  sgt_ncid(:),                &
  sgt_vid(:,:),               &
  sgt_tid(:),                 &
  vel_ncid(:),                &
  vel_vid(:,:),               &
  vel_tid(:)
logical,allocatable,public :: &
  snap_ishere(:),             &
  snap_onfree(:),             &
  sgt_oflag(:),               &
  vel_oflag(:),               &
  sgt_out(:),                 &
  vel_out(:)

character (len=SEIS_STRLEN),public :: pnm_out,pnm_station,  &
  pnm_rest, fnm_rest, fnm_rest_point

!---------------------------------------------
integer pt_ncid,pt_tid,pt_vid(2),pt_sid(1)
integer rest_tinv,run_from_rest
integer :: nt_dyn_rest, nt_dyn_sync, nt_dyn_new

real(SP),dimension(:,:),allocatable :: Vx0,Vz0

!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
subroutine io_init(fnm_input)
  character (len=*),intent(in) :: fnm_input
  integer :: fid

  fid=1001
  open(fid,file=trim(fnm_input),status="old")

  call string_conf(fid,1,'OUTPUT_ROOT',2,pnm_out)
  call string_conf(fid,1,'STATION_ROOT',2,pnm_station)
  !-- log file --
  call string_conf(fid,1,'fnm_log',2,fnm_log)

  !-- restart --
  call string_conf(fid,1,'CHECKPOINT_ROOT',2,pnm_rest)
  call string_conf(fid,1,'checkpoint_tinv',2,rest_tinv)
  call string_conf(fid,1,'run_from_checkpoint',2,run_from_rest)
  call string_conf(fid,1,'urgent_checkpoint',2,fnm_rest_point)
  nt_dyn_rest=0; nt_dyn_sync=0;

  close(fid)

end subroutine io_init

!---------------------------------------------------------------------------

subroutine io_pt_read(fnm_input)
  character (len=*),intent(in) :: fnm_input
  real(SP),dimension(SEIS_GEO) :: xz0,dxz
  integer :: fid,n,m,npt

  fid=1001
  open(fid,file=trim(fnm_input),status="old")
  call string_conf(fid,1,'number_of_recv',2,num_recv)
  call string_conf(fid,1,'number_of_inline',2,num_line)
  call string_conf(fid,1,'tinv_of_seismo',2,pt_tinv)
  call string_conf(fid,1,'topo_hyper_height',2,topo_hyper_height)
  num_pt=num_recv
  do n=1,num_line
    call string_conf(fid,1,trim(io_enum('line_',n)),6,m)
    num_pt=num_pt+m
  end do
  call alloc_pt(num_pt,num_line)
  !-- recv --
  npt=0
  do n=1,num_recv
    npt=npt+1
    do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('recv_',n)),m+1,pt_xz(m,npt))
    end do
    pt_id(:,npt)=(/ n,0 /)
  end do
  !-- inline --
  do n=1,num_line
    do m=1,SEIS_GEO
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1,xz0(m))
      call string_conf(fid,1,trim(io_enum('line_',n)),m+1+SEIS_GEO,dxz(m))
    end do
    call string_conf(fid,1,trim(io_enum('line_',n)),1+1+2*SEIS_GEO,num_line_pt(n))
    do m=1,num_line_pt(n)
      npt=npt+1
      pt_xz(:,npt)=xz0+dxz*(m-1)
      pt_id(:,npt)=(/ m,n /)
    end do
  end do
end subroutine io_pt_read
subroutine alloc_pt(npt,nline)
  integer,intent(in) :: npt,nline
  allocate(pt_xz(SEIS_GEO,npt)); pt_xz=0.0
  allocate(pt_id (       2,npt)); pt_id =0
  allocate(num_line_pt(nline));   num_line_pt=0
end subroutine alloc_pt

subroutine io_pt_import
  integer :: n
  character (len=SEIS_STRLEN) :: filenm
  filenm= get_fnm_station(pnm_station)
  call nfseis_diminfo(trim(filenm),'num_pt',num_pt)
  call nfseis_attget(trim(filenm),'tinv',pt_tinv)
  if (num_pt>0) then
    allocate(pt_indx(SEIS_GEO,num_pt)); pt_indx=0
    call nfseis_varget(trim(filenm),'indx', pt_indx, &
      (/1,1/),(/SEIS_GEO,num_pt/),(/1,1/))
    !do n=1,num_pt
    !   pt_indx(1,n)=inn_i(pt_indx(1,n))
    !   pt_indx(2,n)=inn_j(pt_indx(2,n))
    !   pt_indx(3,n)=inn_k(pt_indx(3,n))
    !end do
  end if
end subroutine io_pt_import

!---------------------------------------------------------------------------

subroutine io_snap_read(fnm_input)
  character (len=*),intent(in) :: fnm_input
  character (len=2) :: varnm
  integer :: fid,n,m

  fid=1001
  open(fid,file=trim(fnm_input),status="old")
  !-- output --
  call string_conf(fid,1,'number_of_snap',2,num_snap)
  call alloc_snap(num_snap)
  do n=1,num_snap
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

subroutine alloc_snap(nsnap)
  integer,intent(in) :: nsnap
  allocate(snap_gsubs(SEIS_GEO,nsnap))
  allocate(snap_gsube(SEIS_GEO,nsnap))
  allocate(snap_gsubt(SEIS_GEO,nsnap))
  allocate(snap_gsubc(SEIS_GEO,nsnap))
  allocate(snap_subs(SEIS_GEO,nsnap))
  allocate(snap_sube(SEIS_GEO,nsnap))
  allocate(snap_subt(SEIS_GEO,nsnap))
  allocate(snap_subc(SEIS_GEO,nsnap))
  allocate(snap_tinv(nsnap))
  allocate(snap_tcnt(nsnap))
  allocate(snap_ishere(nsnap)); snap_ishere=.false.
  allocate(snap_onfree(nsnap)); snap_onfree=.false.
  allocate(vel_oflag(nsnap)); vel_oflag=.false.
  allocate(vel_ncid(nsnap))
  allocate(vel_vid(2,nsnap))
  allocate(vel_tid(nsnap))
  allocate(vel_out(nsnap)); vel_out=.false.
  allocate(sgt_oflag(nsnap)); sgt_oflag=.false.
  allocate(sgt_ncid(nsnap))
  allocate(sgt_vid(2*SEIS_GEO,nsnap))
  allocate(sgt_tid(nsnap))
  allocate(sgt_out(nsnap)); sgt_out=.false.
end subroutine alloc_snap

subroutine io_snap_locate(n_i,n_k)
  integer,intent(in) :: n_i,n_K
  integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
  integer,dimension(SEIS_GEO) :: gsubs,gsube
  integer n

  do n=1,num_snap
    subs=snap_subs(:,n);subc=snap_subc(:,n);subt=snap_subt(:,n)
    sube=snap_sube(:,n)
    ! convert into this thread
    if (      subs(1)<=ngi2 .and. sube(1)>=ngi1 &
      .and. subs(2)<=ngk2 .and. sube(2)>=ngk1 ) then
      snap_ishere(n)=.true.
      if (freenode .and. subs(2)==ngk2) snap_onfree(n)=.true.
      call corr_subse(subs,subc,subt,sube)
      gsubs(1)=out_i(subs(1)); gsubs(2)=out_k(subs(2))
      gsube(1)=out_i(sube(1)); gsube(2)=out_k(sube(2))
      snap_gsubs(:,n)=gsubs
      snap_gsube(:,n)=gsube
      snap_gsubc(:,n)=subc
      snap_subs(1,n)=swmpi_locli(subs(1),n_i)
      snap_subs(2,n)=swmpi_loclk(subs(2),n_k)
      snap_sube(1,n)=swmpi_locli(sube(1),n_i)
      snap_sube(2,n)=swmpi_loclk(sube(2),n_k)
      snap_subc(:,n)=subc
    end if
  end do
  call alloc_V0

  allocate(iovar(nx1:nx2,nz1:nz2)); iovar=0.0

end subroutine io_snap_locate

subroutine alloc_V0
  allocate(Vx0(nx1:nx2,1)); Vx0=0.0
  allocate(Vz0(nx1:nx2,1)); Vz0=0.0
end subroutine alloc_V0
subroutine dealloc_V0
  if (allocated(Vx0)) deallocate(Vx0)
  if (allocated(Vz0)) deallocate(Vz0)
end subroutine dealloc_V0
!---------------------------------------------------------------------------

subroutine corr_indx(subs,subc,subt,sube)
  integer,dimension(SEIS_GEO),intent(inout) :: subs,subc,subt,sube
  ! -1 global first index
  ! 0 index of source center
  ! -2 global index of free surface
  if (subs(1)==-1) then
    subs(1)=ni1
  else
    subs(1)=inn_i(subs(1))
  end if

  if (subs(2)==-2) then
    subs(2)=swmpi_globk(nk2,dims(2)-1)
  elseif (subs(2)==-1) then
    subs(2)=nk1
  else
    subs(2)=inn_k(subs(2))
  end if

  sube=subs+subt*(subc-1)
  if (subc(1)==-1) then
    sube(1)=swmpi_globi(ni2,dims(1)-1)
    subc(1)=(sube(1)-subs(1))/subt(1)+1
  end if
  if (subc(2)==-1) then
    sube(2)=swmpi_globk(nk2,dims(2)-1)
    subc(2)=(sube(2)-subs(2))/subt(2)+1
  end if
end subroutine corr_indx
subroutine corr_subse(subs,subc,subt,sube)
  ! in global mode
  integer,dimension(SEIS_GEO),intent(inout) :: subs,subc,subt,sube
  integer n
  if (ngi1>subs(1)) then
    n=ngi1-subs(1)
    if (mod(n,subt(1))==0) then
      subs(1)=ngi1
    else
      n=(n+subt(1)-1)/subt(1)
      subs(1)=n*subt(1)+subs(1)
    end if
  end if
  if (ngk1>subs(2)) then
    n=ngk1-subs(2)
    if (mod(n,subt(2))==0) then
      subs(2)=ngk1
    else
      n=(n+subt(2)-1)/subt(2)
      subs(2)=n*subt(2)+subs(2)
    end if
  end if
  subc(1)=(ngi2-subs(1))/subt(1)
  subc(2)=(ngk2-subs(2))/subt(2)
  sube(1)=min(sube(1),subc(1)*subt(1)+subs(1))
  sube(2)=min(sube(2),subc(2)*subt(2)+subs(2))
  subc=(sube-subs)/subt+1
end subroutine corr_subse
!---------------------------------------------------------------------------
subroutine io_rest_export(P,Vx,Vz,ntime,is_rest)
  integer,intent(in) :: ntime
  logical,intent(out) :: is_rest
  real(SP),dimension(:,:),intent(in) ::P,Vx,Vz
  integer,dimension(SEIS_GEO) :: subs,subc,subt
  character (len=SEIS_STRLEN) :: filenm
  integer :: fid,n,ierr,ncid,vid(2),sid(1)
  fid=5001
  is_rest=.false.
  if (mod(ntime,10)==0) then
    if (masternode) then
      open(fid,file=trim(fnm_rest_point),status='old',iostat=ierr)
      if (ierr/=0)  &
        call swmpi_except("io_rest_export: error when open "//trim(fnm_rest_point))
      read(fid,*) nt_dyn_rest, nt_dyn_sync,nt_dyn_new
      close(fid)
    end if
    if (nt_dyn_new>ntime) call reset_nt(nt_dyn_new)
  end if

  if (ntime==nt_dyn_rest .or. ntime==nt_dyn_sync .or. mod(ntime,rest_tinv)==0) then
    if (num_pt>0) ierr=nf90_sync(pt_ncid)
    do n=1,num_snap
    if (vel_oflag(n)) ierr=nf90_sync(vel_ncid(n))
    if (sgt_oflag(n)) ierr=nf90_sync(sgt_ncid(n))
    end do
  end if

  if (ntime/=nt_dyn_rest .and. mod(ntime,rest_tinv)/=0) return

  is_rest=.true.

  subs=(/ nx1,nz1 /); subc=(/ nx,nz /); subt=(/ 1,1 /)
  filenm=get_fnm_rest(pnm_rest,ntime)
  call nfseis_grid2d_def(filenm,nx,nz,ncid, &
    "Restart file generated by seis3d_wave" )
  call nfseis_grid2d_defvar(ncid,'Vx' ,vid(1))
  call nfseis_grid2d_defvar(ncid,'Vz' ,vid(2))
  call nfseis_grid2d_defvar(ncid,'P',sid(1))
  call nfseis_grid2d_enddef(ncid)
  call nfseis_put(ncid,vid(1),Vx, subs,subc,subt)
  call nfseis_put(ncid,vid(2),Vz, subs,subc,subt)
  call nfseis_put(ncid,sid(1),P,subs,subc,subt)
  call nfseis_close(ncid)

  if (ntime-rest_tinv>0) then
    filenm=get_fnm_rest(pnm_rest,ntime-rest_tinv)
    call io_delete(filenm)
  end if
end subroutine io_rest_export
subroutine io_rest_import(P,Vx,Vz,ntime)
  integer,intent(out) :: ntime
  real(SP),dimension(:,:),intent(out) :: P,Vx,Vz
  integer,dimension(SEIS_GEO) :: subs,subc,subt
  character (len=SEIS_STRLEN) :: filenm
  if (run_from_rest==0) then
    ntime=0
  else
    ntime=run_from_rest
    subs=(/ nx1,nz1 /); subc=(/ nx,nz /); subt=(/ 1,1 /)
    filenm=get_fnm_rest(pnm_rest,ntime)
    call nfseis_varget( filenm,'Vx',Vx,subs,subc,subt)
    call nfseis_varget( filenm,'Vz',Vz,subs,subc,subt)
    call nfseis_varget( filenm,'P',P,subs,subc,subt)
  end if
end subroutine io_rest_import
!---------------------------------------------------------------------------

subroutine io_seismo_init
  character (len=SEIS_STRLEN) :: filenm
  if (num_pt<=0) return
  filenm=get_fnm_seismo(pnm_out)
  if (run_from_rest==0) then
    call nfseis_seismo_def(filenm,num_pt,pt_ncid,pt_tid,title="seismograms")
    call nfseis_seismo_defvar(pt_ncid,'Vx',pt_vid(1))
    call nfseis_seismo_defvar(pt_ncid,'Vz',pt_vid(2))
    call nfseis_seismo_defvar(pt_ncid,'P',pt_sid(1))
    call nfseis_seismo_enddef(pt_ncid)
  else
    call nfseis_open(filenm,pt_ncid)
    call nfseis_inq_varid(pt_ncid,'time',pt_tid)
    call nfseis_inq_varid(pt_ncid,'Vx',pt_vid(1))
    call nfseis_inq_varid(pt_ncid,'Vz',pt_vid(2))
    call nfseis_inq_varid(pt_ncid,'P',pt_sid(1))
  end if
end subroutine io_seismo_init

subroutine io_seismo_put(Vx,Vz,P,ntime)
  real(SP),dimension(:,:),intent(in) :: Vx,Vz
  real(SP),dimension(:,:),intent(in) :: P
  integer,intent(in) :: ntime
  real(SP) :: t
  integer i,k,n,m,ierr

  if (num_pt<=0) return

  t=ntime*stept
  if (mod(ntime,pt_tinv)==0) then
    m=ntime/pt_tinv
    do n=1,num_pt
      i=pt_indx(1,n);k=pt_indx(2,n)
      call nfseis_put(pt_ncid,pt_tid,t,(/m/),(/1/),(/1/))
      call nfseis_put(pt_ncid,pt_sid(1),P(i,k), &
        (/n,m/),(/1,1/),(/1,1/))
      if(ResultOnColocatedGrid .eqv. .true.) then
        call nfseis_put(pt_ncid,pt_vid(1),(Vx(i,k)+Vx(i-1,k))/2.0, &
          (/n,m/),(/1,1/),(/1,1/))
        call nfseis_put(pt_ncid,pt_vid(2),(Vz(i,k)+Vz(i,k-1))/2.0, &
          (/n,m/),(/1,1/),(/1,1/))
      else
        call nfseis_put(pt_ncid,pt_vid(1),Vx(i,k), &
          (/n,m/),(/1,1/),(/1,1/))
        call nfseis_put(pt_ncid,pt_vid(2),Vz(i,k), &
          (/n,m/),(/1,1/),(/1,1/))
      end if
    end do
  end if

  if (mod(ntime,pt_tinv*100)==0) ierr=nf90_sync(pt_ncid)
end subroutine io_seismo_put

subroutine io_seismo_close
  if (num_pt<=0) return
  call nfseis_close(pt_ncid)
end subroutine io_seismo_close

!---------------------------------------------------------------------------
subroutine io_wave_export(Vx,Vz,ntime,stept)
  real(SP),dimension(:,:),intent(in) :: Vx,Vz
  integer,intent(in) :: ntime
  real(SP),intent(in) :: stept

  integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
  integer,dimension(SEIS_GEO) :: isubs,isube
  integer,dimension(SEIS_GEO+1) :: startput,countput,strideput
  character (len=SEIS_STRLEN) :: filenm
  real(SP) :: t
  integer n,n1

  integer i,k,ierr,si,sk
  real(SP) :: p1(1)

  t=ntime*stept
  !t=cur_time
  do n=1,num_snap
    if (.not. snap_ishere(n)) cycle
    if (.not. vel_out(n)) cycle
    if (mod(ntime,snap_tinv(n))/=0) cycle

    n1=mod(ntime/snap_tinv(n)-1,snap_tcnt(n))
    subs=snap_subs(:,n); subt=snap_subt(:,n); sube=snap_sube(:,n)
    subc=snap_subc(:,n)
    isubs(1) =out_i(subs(1));isubs(2)=out_k(subs(2))
    isube(1) =out_i(sube(1));isube(2)=out_k(sube(2))

    if ( n1==0 ) then
      filenm=get_fnm_snapnode(pnm_out,'vel_',n,ntime,thisid(1),thisid(2))
      call nfseis_snap_def(filenm,vel_ncid(n),vel_tid(n),stept*snap_tinv(n), &
        subc,"snap of velocity feilds")
      call nfseis_snap_attdef(vel_ncid(n), &
        snap_gsubs(:,n), snap_gsubc(:,n), snap_gsubt(:,n), snap_gsube(:,n), &
        isubs, subc, subt, isube)
      call nfseis_snap_defvar(vel_ncid(n),'Vx',vel_vid(1,n))
      call nfseis_snap_defvar(vel_ncid(n),'Vz',vel_vid(2,n))
      call nfseis_snap_enddef(vel_ncid(n))
      vel_oflag(n)=.true.
    elseif ( run_from_rest>0 .and. ntime-run_from_rest<=snap_tinv(n) ) then
      filenm=get_fnm_snapnode(pnm_out,'vel_',n,ntime,thisid(1),thisid(2))
      call nfseis_open(filenm,vel_ncid(n))
      call nfseis_inq_varid(vel_ncid(n),'Vx',vel_vid(1,n))
      call nfseis_inq_varid(vel_ncid(n),'Vz',vel_vid(2,n))
      call nfseis_inq_varid(vel_ncid(n),'time',vel_tid(n))
      vel_oflag(n)=.true.
    end if

    startput=(/1,1,n1+1/);countput=(/subc,1/);strideput=(/1,1,1/)

    call nfseis_put(vel_ncid(n),vel_tid(n),t,(/n1+1/),(/1/),(/1/) )

    if (ResultOnColocatedGrid .eqv. .true.) then

      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(Vx(i,k)+Vx(i-1,k))/2.0
      end do
      end do
      call nfseis_put(vel_ncid(n),vel_vid(1,n), &
        iovar(subs(1):sube(1):subt(1),       &
        subs(2):sube(2):subt(2)),      &
        startput,countput,strideput )
      do k=subs(2),sube(2),subt(2)
      do i=subs(1),sube(1),subt(1)
        iovar(i,k)=(Vz(i,k)+Vz(i,k-1))/2.0
      end do
      end do
      call nfseis_put(vel_ncid(n),vel_vid(2,n),   &
        iovar(subs(1):sube(1):subt(1),         &
        subs(2):sube(2):subt(2)),        &
        startput,countput,strideput )
    else
      call nfseis_put(vel_ncid(n),vel_vid(1,n), &
        Vx(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
      call nfseis_put(vel_ncid(n),vel_vid(2,n), &
        Vz(subs(1):sube(1):subt(1),          &
        subs(2):sube(2):subt(2)),         &
        startput,countput,strideput )
    end if

    if ( n1==snap_tcnt(n)-1 ) then
      call nfseis_close(vel_ncid(n))
      vel_oflag(n)=.false.
    end if

  end do !n
end subroutine io_wave_export

subroutine io_wave_close
  integer :: n
  do n=1,num_snap
    if (vel_oflag(n)) call nfseis_close(vel_ncid(n))
  end do
end subroutine io_wave_close

subroutine io_stress_export(P,ntime,stept)
  real(SP),dimension(:,:),intent(in) :: P
  integer,intent(in) :: ntime
  real(SP),intent(in) :: stept

  integer,dimension(SEIS_GEO) :: subs,subc,subt,sube
  integer,dimension(SEIS_GEO) :: isubs,isube
  integer,dimension(SEIS_GEO+1) :: startput,countput,strideput
  character (len=SEIS_STRLEN) :: filenm
  real(SP) :: t
  integer :: n,n1
  integer i,k,ierr,si,sk
  real(SP) :: p1(1)
  t=(ntime-0.5)*stept
  !t=cur_time
  do n=1,num_snap
    if (.not. snap_ishere(n)) cycle
    if (.not. sgt_out(n)) cycle
    if (mod(ntime,snap_tinv(n))/=0) cycle

    n1=mod(ntime/snap_tinv(n)-1,snap_tcnt(n))
    subs=snap_subs(:,n); subt=snap_subt(:,n); sube=snap_sube(:,n)
    subc=snap_subc(:,n)
    isubs(1) =out_i(subs(1));isubs(2)=out_k(subs(2))
    isube(1) =out_i(sube(1));isube(2)=out_k(sube(2))

    if ( n1==0 ) then
      filenm=get_fnm_snapnode(pnm_out,'sgt_',n,ntime,thisid(1),thisid(2))
      call nfseis_snap_def(filenm,sgt_ncid(n),sgt_tid(n),stept*snap_tinv(n),     &
        subc,"snap of velocity feilds")
      call nfseis_snap_attdef(sgt_ncid(n), &
        snap_gsubs(:,n), snap_gsubc(:,n), snap_gsubt(:,n), snap_gsube(:,n), &
        isubs, subc, subt, isube)
      call nfseis_snap_defvar(sgt_ncid(n),'P',sgt_vid(1,n))
      call nfseis_snap_enddef(sgt_ncid(n))
      sgt_oflag(n)=.true.
    elseif ( run_from_rest>0 .and. ntime-run_from_rest<=snap_tinv(n) ) then
      filenm=get_fnm_snapnode(pnm_out,'sgt_',n,ntime,thisid(1),thisid(2))
      call nfseis_open(filenm,sgt_ncid(n))
      call nfseis_inq_varid(sgt_ncid(n),'P',sgt_vid(1,n))
      call nfseis_inq_varid(sgt_ncid(n),'time',sgt_tid(n))
      sgt_oflag(n)=.true.
    end if

    startput=(/1,1,n1+1/);countput=(/subc,1/);strideput=(/1,1,1/)

    call nfseis_put(sgt_ncid(n),sgt_tid(n),t,(/n1+1/),(/1/),(/1/) )

    call nfseis_put(sgt_ncid(n),sgt_vid(1,n), &
       P(subs(1):sube(1):subt(1),         &
      subs(2):sube(2):subt(2)),        &
      startput,countput,strideput )
  end do

end subroutine io_stress_export
subroutine io_stress_close
  integer :: n
  do n=1,num_snap
    if (sgt_oflag(n)) call nfseis_close(sgt_ncid(n))
  end do
end subroutine io_stress_close

!---------------------------------------------------------------------------
function get_fnm_station(pnm_sta) result(filenm)
  character (len=*),intent(in) :: pnm_sta
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm_sta)                    &
    //'/'//'station'                   &
    //'.nc'
end function get_fnm_station
function get_fnm_seismo(pnm) result(filenm)
  character (len=*),intent(in) :: pnm
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm    )                          &
    //'/'//'seismo'                          &
    //'.nc'
end function get_fnm_seismo
function get_fnm_snapnode(pnm,prefix,n,ntime,n_i,n_k) result(filenm)
  integer,intent(in) :: n,ntime,n_i,n_k
  character (len=*),intent(in) :: pnm,prefix
  character (len=SEIS_STRLEN) :: filenm
  integer n0
  n0=(ntime+snap_tinv(n)*snap_tcnt(n)-1)/(snap_tinv(n)*snap_tcnt(n))
  filenm=trim(pnm)                              &
    //'/'//trim(io_enum(prefix,n))           &
    //'_n'//trim(io_out_pattern(n0,5))       &
    //'.nc'
end function get_fnm_snapnode
function get_fnm_snapnode_n(pnm,prefix,id,n0,n_i,n_k) result(filenm)
  character (len=*),intent(in) :: pnm,prefix
  integer,intent(in) :: id,n0,n_i,n_k
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm)                              &
    //'/'//trim(io_enum(prefix,id))          &
    //'_n'//trim(io_out_pattern(n0,5))       &
    //'.nc'
end function get_fnm_snapnode_n
function get_fnm_rest(pnm,ntime) result(filenm)
  integer,intent(in) :: ntime
  character (len=*),intent(in) :: pnm
  character (len=SEIS_STRLEN) :: filenm
  filenm=trim(pnm)//'/'                          &
    //'rest_t'//trim(io_out_pattern(ntime,5)) &
    //'.nc'
end function get_fnm_rest

function io_enum(prefix,num) result(ioname)
  character (len=*),intent(in) :: prefix
  integer,intent(in) :: num
  character (len=SEIS_STRLEN) :: ioname
  character (len=SEIS_STRLEN) :: str
  write(str,"(i3.3)") num
  ioname=trim(prefix)//trim(str)
end function io_enum

function io_out_pattern(num,width) result(ioname)
  integer,intent(in) :: num
  integer,optional,intent(in) :: width
  character (len=SEIS_STRLEN) :: ioname,str,fmt_str
  if (present(width)) then
    write(str,"(i5)") width
    fmt_str="(i"//trim(str)//"."//trim(str)//")"
  else
    fmt_str="(i4.4)"
  end if
  write(ioname,fmt_str) num
end function io_out_pattern

subroutine io_delete(filenm)
  character (len=*),intent(in) :: filenm
  integer fid
  fid=5001
  open(fid,file=trim(filenm),status='unknown')
  close(fid,status='delete')
end subroutine io_delete

end module io_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
