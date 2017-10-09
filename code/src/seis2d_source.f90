

program seis2d_source

!----------------------------------------------------------------------------
! This program generate sources.
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------






































































use netcdf
use constants_mod
use string_mod
use math_mod
use nfseis_mod
use para_mod
use grid_mod
use media_mod
use src_mod
use nompi_mod
implicit none

real(SP) :: src_hyper_height

logical,allocatable :: force_flag(:),moment_flag(:),resrc_flag(:)
real(SP),allocatable :: force_axis(:,:), moment_axis(:,:), resrc_axis(:,:)
character (len=SEIS_STRLEN) :: filenm
integer :: p(1)
real(SP) :: xmin,xmax,zmin,zmax
real(SP) :: x0,z0
integer :: i,k,n_i,n_k
integer :: n,m,npt,nfrc,nmom,gi,gk,si,sk,nresrc,nresrc_point
logical :: iflag
integer :: ncid,indxid,axisid,siftid,t0id,stftid,stffid
integer :: fxid,fzid,mppid,vxid,vzid,resrc_pointid

!----------------------------------------------------------------------!

call get_conf_name(fnm_conf)! !fnm_conf='SeisFD3D.conf'

call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call grid_fnm_init(fnm_conf)
call media_fnm_init(fnm_conf)
call src_fnm_init(fnm_conf)

! load grid
call grid_coord_alloc

! source
call read_src_para(fnm_src_conf,fnm_conf)

! find src index
n_i=0; n_k=0;
write(*,"(i10,i2,a,2(i2))") n_i,n_k, ' of ',dims
call swmpi_set_gindx(n_i,n_k)
call grid_coord_import

xmin=minval(x); xmax=maxval(x)
zmin=minval(z); zmax=maxval(z)
print *, 'xmin=',xmin,'xmax=',xmax
print *, 'zmin=',zmin,'zmax=',zmax

!force
do npt=1,num_force
  if (force_flag(npt)) cycle

  x0=force_axis(1,npt);z0=force_axis(2,npt)
  if (x0<xmin .or. x0>xmax ) cycle

  p=minloc(abs(x-x0)); i=loct_i(p(1))

  if( z0 >=zmin .and. z0 <=zmax ) then
    p=minloc(abs(z-z0)); k=loct_k(p(1))
  elseif (z0>src_hyper_height .and. n_k==dims(2)-1) then
    k=nk2; z0=z(k)
  else
    cycle
  end if

  if ( i<ni1 .or. i>ni2 .or. k<nk1 .or. k>nk2 ) cycle

  write(*,"(i10,a,i10,1(i2),a,2(i2))") npt, ' in', n_i,n_k, ' of ',dims
  force_flag(npt)=.true.
  gi=swmpi_globi(i,n_i); gk=swmpi_globk(k,n_k)
  force_indx(:,npt)=(/ gi,gk /)
  force_axis(2,npt)=z0

  ! calculate shift
  force_shift(:,npt)=0.0_SP
  if (x0>x(i)) then
    force_shift(1,npt)=(x0-x(i))/(x(i+1)-x(i))
  else
    force_shift(1,npt)=(x0-x(i))/(x(i)-x(i-1))
  end if
  if (z0>z(k)) then
    force_shift(2,npt)=(z0-z(k))/(z(k+1)-z(k))
  else
    force_shift(2,npt)=(z0-z(k))/(z(k)-z(k-1))
  end if
end do !n


!moment
do npt=1,num_moment
  if (moment_flag(npt)) cycle

  x0=moment_axis(1,npt);z0=moment_axis(2,npt)
  if (x0<xmin .or. x0>xmax) cycle

  p=minloc(abs(x-x0)); i=loct_i(p(1))

  if( z0 >=zmin .and. z0 <=zmax ) then
    p=minloc(abs(z-z0)); k=loct_k(p(1))
  elseif (z0>src_hyper_height .and. n_k==dims(2)-1) then
    k=nk2; z0=z(k)
  else
    cycle
  end if

  if ( i<ni1 .or. i>ni2 .or. k<nk1 .or. k>nk2 ) cycle

  write(*,"(i10,a,i10,1(i2),a,2(i2))") npt, ' in', n_i,n_k, ' of ',dims
  moment_flag(npt)=.true.
  gi=swmpi_globi(i,n_i); gk=swmpi_globk(k,n_k)
  moment_indx(:,npt)=(/ gi,gk /)
  moment_axis(2,npt)=z0

  ! calculate shift
  moment_shift(:,npt)=0.0_SP
  if (x0>x(i)) then
    moment_shift(1,npt)=(x0-x(i))/(x(i+1)-x(i))
  else
    moment_shift(1,npt)=(x0-x(i))/(x(i)-x(i-1))
  end if
  if (z0>z(k)) then
    moment_shift(2,npt)=(z0-z(k))/(z(k+1)-z(k))
  else
    moment_shift(2,npt)=(z0-z(k))/(z(k)-z(k-1))
  end if
end do !n
!resrc
do npt=1,num_resrc
    if (resrc_flag(npt)) cycle
    x0 = resrc_axis(1,npt); z0 = resrc_axis(2,npt); 
!    print *,'x0=',x0
!    print *,'z0=',z0
    if (x0<xmin .or. x0>xmax) cycle
    p=minloc(abs(x-x0)); i=loct_i(p(1))
    if( z0 >=zmin .and. z0 <=zmax ) then
       p=minloc(abs(z-z0)); k=loct_k(p(1))
    elseif (z0>src_hyper_height .and. n_k==dims(2)-1) then
       k=nk2; z0=z(k)
    else
       cycle
    end if
     if ( i<ni1 .or. i>ni2 .or. k<nk1 .or. k>nk2 ) cycle
    write(*,"(i10,a,i10,1(i2),a,2(i2))") npt, ' in', n_i,n_k, ' of',dims
    resrc_flag(npt)=.true.
    gi=swmpi_globi(i,n_i);gk=swmpi_globk(k,n_k)
    resrc_indx(:,npt)=(/ gi,gk /)
    resrc_axis(2,npt)=z0
end do

! check
iflag=.false.
do n=1,num_force
  if ( .not. force_flag(n)) then
    iflag=.true.
    print *, 'n=',n,'loct=',force_axis(:,n),'indx=',force_indx(:,n)
  end if
end do
do n=1,num_moment
  if ( .not. moment_flag(n)) then
    iflag=.true.
    print *, 'n=',n,'loct=',moment_axis(:,n),'indx=',moment_indx(:,n)
  end if
end do
do n=1,num_resrc
   if ( .not. resrc_flag(n)) then
      iflag=.true.
      print *, 'n=',n,'loct=',resrc_axis(:,n),'indx=',resrc_indx(:,n)
   end if
end do
if (iflag) then
  print *, 'there are some source points out of the computational domain'
  stop 1
end if

nfrc=0
do n=1,num_force
  gi=force_indx(1,n);gk=force_indx(2,n)
  if (      gi>=ngx1 .and. gi<=ngx2                          &
    .and. gk>=ngz1 .and. gk<=ngz2 ) nfrc=nfrc+1
end do
nmom=0
do n=1,num_moment
  gi=moment_indx(1,n);gk=moment_indx(2,n)
  if (      gi>=ngx1 .and. gi<=ngx2                          &
    .and. gk>=ngz1 .and. gk<=ngz2 ) nmom=nmom+1
end do
nresrc=0
do n=1,num_resrc
  gi=resrc_indx(1,n);gk=resrc_indx(2,n)
  if (      gi>=ngx1 .and. gi<=ngx2                          &
    .and. gk>=ngz1 .and. gk<=ngz2 ) nresrc=nresrc+1
end do

filenm=src_fnm_get()
call srcnode_skel(filenm,nfrc,ntwin_force,nmom,nresrc,ntwin_moment)
call nfseis_open(filenm,ncid)

!force
if (nfrc>0) then
  call nfseis_inq_varid(ncid,'force_indx',indxid)
  call nfseis_inq_varid(ncid,'force_axis',axisid)
  call nfseis_inq_varid(ncid,'force_shift',siftid)
  call nfseis_inq_varid(ncid,'force_start_time',t0id)
  call nfseis_inq_varid(ncid,'force_stf_time',stftid)
  call nfseis_inq_varid(ncid,'force_stf_freq',stffid)
  call nfseis_inq_varid(ncid,'Fx',fxid)
  call nfseis_inq_varid(ncid,'Fz',fzid)
  call nfseis_put(ncid,stftid,frcstf_time,(/1/),(/ntwin_force/),(/1/))
  call nfseis_put(ncid,stffid,frcstf_freq,(/1/),(/ntwin_force/),(/1/))
  m=0
  do n=1,num_force
    gi=force_indx(1,n);gk=force_indx(2,n)
    if (      gi>=ngx1 .and. gi<=ngx2                          &
      .and. gk>=ngz1 .and. gk<=ngz2 ) then
      m=m+1
      call nfseis_put(ncid,fxid,ForceX(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
      call nfseis_put(ncid,fzid,ForceZ(:,n),(/1,m/),(/ntwin_force,1/),(/1,1/))
      call nfseis_put(ncid,axisid,force_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_put(ncid,siftid,force_shift(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_put(ncid,t0id,force_t0(n),(/m/),(/1/),(/1/))
      call nfseis_put(ncid,indxid,         &
        !-- (/ swmpi_locli(gi,n_i),         &
      !--    swmpi_loclk(gk,n_k) /),      &
      (/ gi, gk /),      &
        (/1,m/),(/SEIS_GEO,1/),(/1,1/))
    end if
  end do
end if
!moment
if (nmom>0) then
  call nfseis_inq_varid(ncid,'moment_indx',indxid)
  call nfseis_inq_varid(ncid,'moment_axis',axisid)
  call nfseis_inq_varid(ncid,'moment_shift',siftid)
  call nfseis_inq_varid(ncid,'moment_start_time',t0id)
  call nfseis_inq_varid(ncid,'moment_rate_time',stftid)
  call nfseis_inq_varid(ncid,'moment_rate_freq',stffid)
  call nfseis_inq_varid(ncid,'PP',mppid)
  call nfseis_put(ncid,stftid,momstf_time,(/1/),(/ntwin_moment/),(/1/))
  call nfseis_put(ncid,stffid,momstf_freq,(/1/),(/ntwin_moment/),(/1/))
  m=0
  do n=1,num_moment
    gi=moment_indx(1,n);gk=moment_indx(2,n)
    if (      gi>=ngx1 .and. gi<=ngx2                          &
      .and. gk>=ngz1 .and. gk<=ngz2 ) then
      m=m+1
      call nfseis_put(ncid,mppid,MomPP(:,n),(/1,m/),(/ntwin_moment,1/),(/1,1/))
      call nfseis_put(ncid,axisid,moment_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_put(ncid,siftid,moment_shift(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_put(ncid,t0id,moment_t0(n),(/m/),(/1/),(/1/))
      call nfseis_put(ncid,indxid,         &
        !-- (/ swmpi_locli(gi,n_i),         &
      !--    swmpi_loclk(gk,n_k) /),      &
      (/ gi, gk /),      &
        (/1,m/),(/SEIS_GEO,1/),(/1,1/))
    end if
  end do
end if
!resrc 
if (nresrc>0) then
     call nfseis_inq_varid(ncid,'resrc_indx',indxid)
     call nfseis_inq_varid(ncid,'resrc_axis',axisid)
     call nfseis_inq_varid(ncid,'Vx',vxid)
     call nfseis_inq_varid(ncid,'Vz',vzid)
     m=0
     do n=1,num_resrc
        gi=resrc_indx(1,n);gk=resrc_indx(2,n)
        if(    gi>=ngx1  .and.  gi<=ngx2                           &
           .and. gk>=ngz1 .and. gk<=ngz2 ) then
       m = m+1
      call nfseis_put(ncid,axisid,resrc_axis(:,n),(/1,m/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_put(ncid,vxid,resrc_Vx(:,n),(/1,m/),(/ntt,1/),(/1,1/))
      call nfseis_put(ncid,vzid,resrc_Vz(:,n),(/1,m/),(/ntt,1/),(/1,1/))
      call nfseis_put(ncid,indxid,          &
            (/ gi, gk /),      &
            (/1,m/),(/SEIS_GEO,1/),(/1,1/))
       end if
    end do
end if

call nfseis_close(ncid)


call src_destroy
call dealloc_local

!----------------------------------------------------------------------!
contains
!----------------------------------------------------------------------!

!*************************************************************************
!*                 PART-I  src alloc and dealloc                         *
!*************************************************************************

subroutine alloc_force(npt)
  integer,intent(in) :: npt
  allocate(force_flag(npt)); force_flag=.false.
  allocate(force_axis(SEIS_GEO,npt)); force_axis=0.0_SP
end subroutine alloc_force
subroutine alloc_moment(npt)
  integer,intent(in) :: npt
  allocate(moment_flag(npt)); moment_flag=.false.
  allocate(moment_axis(SEIS_GEO,npt)); moment_axis=0.0_SP
end subroutine alloc_moment
subroutine alloc_resrc(npt)
  integer,intent(in) :: npt
  allocate(resrc_flag(npt)); resrc_flag=.false.
  allocate(resrc_axis(SEIS_GEO,npt)); resrc_axis=0.0_SP
end subroutine alloc_resrc

subroutine dealloc_local
  if (allocated(force_flag)) deallocate(force_flag)
  if (allocated(moment_flag)) deallocate(moment_flag)
  if (allocated(resrc_flag)) deallocate(resrc_flag)
end subroutine dealloc_local

!*************************************************************************
!*                    PART-II  src init and distrib                      *
!*************************************************************************

subroutine read_src_para(fnm_src_conf,fnm_conf)
  character (len=*),intent(in) :: fnm_conf,fnm_src_conf
  character (len=SEIS_STRLEN) :: mommech,stf_type,str
  character (len=SEIS_STRLEN) :: pnm_src,fnm
  real(DP) :: strike,dip,rake
  real(SP),dimension(:,:),allocatable :: A
  real(SP) :: f0,m0
  integer fid1,k1,k2
  integer fid,n,m

  fid=1001
  open(fid,file=trim(fnm_src_conf),status="old")

  call string_conf(fid,1,'src_hyper_height',2,src_hyper_height)

  !force
  call string_conf(fid,1,"number_of_force_source",2,num_force)
  if (num_force>=1) then
    !call string_conf(fid,1,"force_stf_outer",2,force_stf_outer)
    !if (force_stf_outer==1) then
    !  call string_conf(fid,1,"force_outer_fnm",2,force_outer_fnm)
    !end if
    call string_conf(fid,1,"force_stf_window",2,ntwin_force)
    call string_conf(fid,1,"force_stf_type",2,stf_type)
    frcstf_id=stf_name2id(trim(stf_type))
    call src_alloc_force(num_force,ntwin_force)
    call alloc_force(num_force)
    do m=1,ntwin_force
      call string_conf(fid,1,"force_stf_timefactor",m+1,frcstf_time(m))
      call string_conf(fid,1,"force_stf_freqfactor",m+1,frcstf_freq(m))
    end do
    call string_conf(fid,1,"<anchor_force>",1,str)
    do n=1,num_force
    do m=1,ntwin_force
      read(fid,*) force_axis(:,n),force_t0(n), &
        f0,ForceX(m,n),ForceZ(m,n)
      ForceX(m,n)=ForceX(m,n)*f0
      ForceZ(m,n)=ForceZ(m,n)*f0
    end do
    end do
  end if
  !moment
  call string_conf(fid,1,"number_of_moment_source",2,num_moment)
  if (num_moment>=1) then
    !call string_conf(fid,1,"moment_rate_outer",2,moment_rate_outer)
    !if (moment_rate_outer==1) then
    !  call string_conf(fid,1,"moment_outer_fnm",2,moment_outer_fnm)
    !end if
    call string_conf(fid,1,"moment_rate_window",2,ntwin_moment)
    call string_conf(fid,1,"moment_rate_type",2,stf_type)
    momstf_id=stf_name2id(trim(stf_type))
    call src_alloc_moment(num_moment,ntwin_moment)
    call alloc_moment(num_moment)
    do m=1,ntwin_moment
      call string_conf(fid,1,"moment_rate_timefactor",m+1,momstf_time(m))
      call string_conf(fid,1,"moment_rate_freqfactor",m+1,momstf_freq(m))
    end do
    call string_conf(fid,1,"moment_mech_input",2,mommech)
    call string_conf(fid,1,"<anchor_moment>",1,str)
    if (trim(mommech)=='moment') then
      do n=1,num_moment
      do m=1,ntwin_moment
        read(fid,*) moment_axis(:,n),moment_t0(n),m0, &
          MomPP(m,n)
        MomPP(m,n)=m0*MomPP(m,n)
      end do
      end do
    else
      do n=1,num_moment
      do m=1,ntwin_moment
     !  'rate'不可以用
     !   read(fid,*) moment_axis(:,n),moment_t0(n),m0,strike,dip,rake
     !   call angle2moment(strike,dip,rake, &
     !     MomTxx(m,n),MomTzz(m,n),MomTxz(m,n))
     !   MomTxx(m,n)=m0*MomTxx(m,n)
     !   MomTzz(m,n)=m0*MomTzz(m,n)
     !   MomTxz(m,n)=m0*MomTxz(m,n)
      end do
      end do
    end if
  end if
  !resrc
  call string_conf(fid,1,"number_of_resrc",2,num_resrc)
  call string_conf(fid,1,'time',2,ntt)
  if (num_resrc>=1) then
      call src_alloc_resrc(ntt,num_resrc)
      call alloc_resrc(num_resrc)
      call string_conf(fid,1,"<anchor_resrc>",1,str)
      k=1
  do n=1,num_resrc
      read(fid,*) resrc_axis(:,n)
      x0=resrc_axis(1,n); z0=resrc_axis(2,n)
        !call swmpi_change_fnm(n_i,n_k)
        !call swmpi_set_gindx(n_i,n_k)
        !call grid_coord_import(n_i,n_k)
             xmin=minval(x); xmax=maxval(x)
             zmin=minval(z); zmax=maxval(z)
             !print *, 'x0=', x0
             !print *, 'z0=', z0
             !print *, 'xmax=', xmax
             !print *, 'zmax=', zmax
    !     k = A(k1,k2)
         fid1 = 1002
         open(fid1,file=trim(fnm_conf),status='old')
       call string_conf(fid1,1,'SOURCE_ROOT',2,pnm_src)       ! pnm_src--SOURCE_ROOT = ./input
       fnm=trim(pnm_src)//'/'//'seismo'//'_reverse'//'.nc'
       call nfseis_varget(fnm,'Vx_reverse',resrc_Vx(:,n),(/k,1/),(/1,ntt/),(/1,1/))
       call nfseis_varget(fnm,'Vz_reverse',resrc_Vz(:,n),(/k,1/),(/1,ntt/),(/1,1/))
       k=k+1
    !   A(k1,k2) = k
       close(fid1)
  end do
  end if

  close(fid)
end subroutine read_src_para

subroutine angle2moment(strike,dip,rake,Mxx,Mzz,Mxz)
  real(DP),intent(in) :: strike,dip,rake
  real(SP),intent(out) :: Mxx,Mzz,Mxz

  real(DP) :: strike_pi,dip_pi,rake_pi
  real(SP) :: M11,M22,M33,M12,M13,M23

  dip_pi=dip/180.0_DP*PI
  strike_pi=strike/180.0_DP*PI
  rake_pi=rake/180.0_DP*PI
  !in Aki and Richard's
  M11=-(sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)   &
    +sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi)**2)
  M22=sin(dip_pi)*cos(rake_pi)*sin(2.0_DP*strike_pi)     &
    -sin(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi)**2
  !Mzz=sin(2.0*dip_pi)*sin(rake_pi)
  M33=-(M11+M22)
  M12=sin(dip_pi)*cos(rake_pi)*cos(2.0_DP*strike_pi)     &
    +0.5_DP*sin(2.0_DP*dip_pi)*sin(rake_pi)*sin(2.0_DP*strike_pi)
  M13=-(cos(dip_pi)*cos(rake_pi)*cos(strike_pi)       &
    +cos(2.0_DP*dip_pi)*sin(rake_pi)*sin(strike_pi))
  M23=-(cos(dip_pi)*cos(rake_pi)*sin(strike_pi)       &
    -cos(2.0_DP*dip_pi)*sin(rake_pi)*cos(strike_pi))
  !Mxz=-Mxz;Mxy=-Mxy !for upward positive z axis
  Mxx= M22; Mzz=M33; Mxz=-M23
end subroutine angle2moment

subroutine srcnode_skel(filenm,nfrc,ntwfrc,nmom,nresrc,ntwmom)
  character (len=*),intent(in) :: filenm
  integer,intent(in) :: nfrc,ntwfrc,nmom,nresrc,ntwmom

  integer :: ncid,ierr,oldMode
  integer :: geoid,vid,nmomid,ntwmomid,nfrcid,ntwfrcid,nresrcid,nttid

  ierr=nf90_create( path=trim(filenm), cmode= nf90_clobber, ncid = ncid )
  call nfseis_except(ierr,'srcnode_init:'//trim(filenm))
  ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
  call nfseis_except(ierr,'set_fill in srcnode_skel')
  ! -- define dim
  ierr=nf90_def_dim(ncid,'geo_dimension',SEIS_GEO,geoid)
  call nfseis_except(ierr,'geodim dim in srcnode_skel')
  ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_force",nfrc)
  !ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_outer",force_stf_outer)
  ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_moment",nmom)
  !ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_rate_outer",moment_rate_outer)
  ierr=nf90_put_att(ncid,NF90_GLOBAL,"number_of_resrc",nresrc)
  !force
  if (nfrc>0) then
    ierr=nf90_def_dim(ncid,'force_point',nfrc,nfrcid)
    call nfseis_except(ierr,'force_point dim in srcnode_skel')
    ierr=nf90_def_dim(ncid,'force_time_window' ,ntwfrc,ntwfrcid)
    call nfseis_except(ierr,'force_time_window dim in srcnode_skel')
    ! -- define variable
    ierr=nf90_def_var(ncid,'Fx', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
    ierr=nf90_def_var(ncid,'Fz', SEISNC_DATATYPE,(/ ntwfrcid, nfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_indx', nf90_int,(/ geoid, nfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_axis', SEISNC_DATATYPE,(/ geoid, nfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_shift', SEISNC_DATATYPE,(/ geoid, nfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_start_time', SEISNC_DATATYPE,(/ nfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_stf_time', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
    ierr=nf90_def_var(ncid,'force_stf_freq', SEISNC_DATATYPE,(/ ntwfrcid /),vid )
    ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_id",frcstf_id)
    ierr=nf90_put_att(ncid,NF90_GLOBAL,"force_stf_type", &
      trim(stf_id2name(frcstf_id)) )
  end if
  !moment
  if (nmom>0) then
    ierr=nf90_def_dim(ncid,'moment_point',nmom,nmomid)
    call nfseis_except(ierr,'moment_point dim in srcnode_skel')
    ierr=nf90_def_dim(ncid,'moment_time_window' ,ntwmom,ntwmomid)
    call nfseis_except(ierr,'moment_time_window dim in srcnode_skel')
    ! -- define variable
    ierr=nf90_def_var(ncid,'PP',SEISNC_DATATYPE,(/ntwmomid,nmomid/),vid)
    ierr=nf90_def_var(ncid,'moment_indx', nf90_int,(/ geoid, nmomid /),vid )
    ierr=nf90_def_var(ncid,'moment_axis', SEISNC_DATATYPE,(/ geoid, nmomid /),vid )
    ierr=nf90_def_var(ncid,'moment_shift', SEISNC_DATATYPE,(/ geoid, nmomid /),vid )
    ierr=nf90_def_var(ncid,'moment_start_time', SEISNC_DATATYPE,(/ nmomid /),vid )
    ierr=nf90_def_var(ncid,'moment_rate_time', SEISNC_DATATYPE,(/ ntwmomid /),vid )
    ierr=nf90_def_var(ncid,'moment_rate_freq', SEISNC_DATATYPE,(/ ntwmomid /),vid )
    ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_rate_id",momstf_id)
    ierr=nf90_put_att(ncid,NF90_GLOBAL,"moment_rate_type", &
      trim(stf_id2name(momstf_id)) )
  end if
  !resrc
  if (nresrc>0) then
    ierr=nf90_def_dim(ncid,'time',ntt,nttid)
    call nfseis_except(ierr,'time in srcnode_skel')
    ierr=nf90_def_dim(ncid,'resrc_point' ,nresrc,nresrcid)
    call nfseis_except(ierr,'resrc_point dim in srcnode_skel')
    ! -- define variable
    ierr=nf90_def_var(ncid,'Vx',SEISNC_DATATYPE,(/nttid,nresrcid/),vid)
    ierr=nf90_def_var(ncid,'Vz',SEISNC_DATATYPE,(/nttid,nresrcid/),vid)
    ierr=nf90_def_var(ncid,'resrc_indx', nf90_int,(/ geoid, nresrcid /),vid )
    ierr=nf90_def_var(ncid,'resrc_axis', SEISNC_DATATYPE,(/ geoid, nresrcid /),vid )
  end if
  ierr=nf90_enddef(ncid)
  call nfseis_except(ierr,'enddef in srcnode_skel')
  ierr=nf90_close(ncid)
  call nfseis_except(ierr,'file close in srcnode_skel')
end subroutine srcnode_skel

end program seis2d_source

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
