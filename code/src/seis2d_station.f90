

program seis2d_station

!----------------------------------------------------------------------------
! This program generate stations
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------


use constants_mod
use string_mod, only : string_conf
use para_mod
use math_mod
use nompi_mod
use nfseis_mod
use grid_mod
use io_mod

implicit none

!-----------------------------------------------------------------------------

call get_conf_name(fnm_conf)

call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call swmpi_set_gindx(0,0)

call grid_fnm_init(fnm_conf)
call grid_coord_alloc

call io_init(fnm_conf)
call io_pt_read(fnm_conf)

call grid_coord_import
call receiver_locate

call grid_destroy

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

subroutine receiver_locate
  character (len=SEIS_STRLEN) :: filenm
  integer,dimension(1) :: p
  integer npt,n,i,k,gi,gk,n_k
  real(SP) :: x0,z0
  npt=0
  n_k=0
  filenm=get_fnm_station(pnm_station)
  call station_skel(filenm,pt_tinv)
  do n=1,num_pt
    x0=pt_xz(1,n);z0=pt_xz(2,n)
    if (z0<=topo_hyper_height) then
      p=minloc(abs(x-x0)); i=p(1)
      p=minloc(abs(z-z0)); k=p(1)
    elseif (n_k==dims(2)-1) then
      p=minloc(abs(x-x0)); i=p(1)
      k=nk2
    else
      i=nx1;k=nz1
    end if
    if ( i>=ni1 .and. i<=ni2      &
      .and. k>=nk1 .and. k<=nk2 ) then
      npt=npt+1
      !-- gi=swmpi_globi(i,n_i)
      !-- gk=swmpi_globk(k,n_k)
      gi=i; gk=k
      call nfseis_varput(filenm,'indx',(/i,k/),                          &
        (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'gindx',(/out_i(gi),out_k(gk)/), &
        (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'coord',(/x0,z0/),                      &
        (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'grid',(/x(i),z(k)/),                 &
        (/1,npt/),(/SEIS_GEO,1/),(/1,1/))
      call nfseis_varput(filenm,'id',pt_id(:,n),           &
        (/1,npt/),(/2,1/),(/1,1/))
    end if
  end do
end subroutine receiver_locate

subroutine station_skel( filenm,tinv,title )
  character (len=*),intent(in) :: filenm
  integer,intent(in) :: tinv
  character (len=*),optional,intent(in) :: title
  integer :: ierr,ncid,gdimid,ndimid,oldMode
  integer :: indxid,gindxid,coordid,loctid,idid,twoelemid
  !--
  ierr=nf90_create( path = trim(filenm), cmode= nf90_clobber, ncid = ncid )
  call nfseis_except(ierr,'station_skel:'//trim(filenm))
  !--
  ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
  call nfseis_except(ierr,'set_fill in station_skel')
  !-- define dim
  ierr=nf90_def_dim(ncid, "num_pt", nf90_unlimited, ndimid)
  call nfseis_except(ierr,'num_pt dim in station_skel')
  ierr=nf90_def_dim(ncid, "geo_dim", SEIS_GEO, gdimid)
  call nfseis_except(ierr,'geo_dim dim in station_skel')
  ierr=nf90_def_dim(ncid, "twoelem", 2, twoelemid)
  call nfseis_except(ierr,'twoelem dim in station_skel')
  ! -- define variable
  ierr=nf90_def_var(ncid, 'indx', nf90_int, (/ gdimid, ndimid /), indxid )
  call nfseis_except(ierr,'indx var in station_skel')
  ierr=nf90_def_var(ncid, 'gindx', nf90_int, (/ gdimid, ndimid /), gindxid )
  call nfseis_except(ierr,'gindx var in station_skel')
  ierr=nf90_def_var(ncid, 'coord', SEISNC_DATATYPE, (/ gdimid,ndimid /), coordid )
  call nfseis_except(ierr,'coord var in station_skel')
  ierr=nf90_def_var(ncid, 'grid', SEISNC_DATATYPE, (/ gdimid,ndimid /), loctid )
  call nfseis_except(ierr,'grid var in station_skel')
  ierr=nf90_def_var(ncid, 'id', nf90_int, (/ twoelemid,ndimid /), idid )
  call nfseis_except(ierr,'id var in station_skel')
  !-- define global attribute
  ierr=nf90_put_att(ncid,NF90_GLOBAL,"tinv",tinv )
  call nfseis_except(ierr,'tinv att in station_skel')
  if (present(title) ) then
    ierr=nf90_put_att(ncid,NF90_GLOBAL,"title",title )
    call nfseis_except(ierr,'title att in station_skel')
  end if
  !--
  ierr=nf90_enddef(ncid)
  call nfseis_except(ierr,'enddef in station_skel')
  !-- 
  ierr=nf90_close(ncid)
  call nfseis_except(ierr,'file close in station_skel')
end subroutine station_skel

end program seis2d_station

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
