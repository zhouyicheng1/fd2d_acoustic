

program seis2d_grid

!----------------------------------------------------------------------------
! This program generates grid coordinate.
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!----------------------------------------------------------------------------






































































use constants_mod
use string_mod
use para_mod
use math_mod
use nompi_mod
use grid_mod
use nfseis_mod

implicit none

!-----------------------------------------------------------------------------
type STRUCT_TOPO
  real(SP),dimension(:),pointer :: x
  real(SP),dimension(:),pointer :: z
end type STRUCT_TOPO

type STRUCT_GRID
  integer :: ni,nk
  real(SP),dimension(:),pointer :: x,z
end type STRUCT_GRID

type (STRUCT_TOPO) :: P
type (STRUCT_GRID) :: G

!free surface topography
real(SP),allocatable :: topogrid(:),gtopoz(:)
real(SP),allocatable :: gx(:),gy(:),gz(:),gzhalf(:)
integer,allocatable :: gfindx(:)
logical,allocatable :: gbit(:),gbit_mask(:)

character (len=SEIS_STRLEN) :: filenm
integer i,k,n,n_i,n_k

!-----------------------------------------------------------------------------

!fnm_grid_conf='SeisGrid.conf'

call get_conf_name(fnm_conf)

call para_init(fnm_conf)
call swmpi_init(fnm_conf)
call swmpi_set_gindx(0,0)

call grid_fnm_init(fnm_conf)
call grid_coord_alloc

!-----------------------------------------------------------------------------
call init_grid(fnm_grid_conf)

!-----------------------------------------------------------------------------
call alloc_local

gx(ni1:NTPI+2)=G%x;gz(nk1:NTPK+2)=G%z
do n=1,2
gx(ni1-n)=gx(ni1-n+1)+gx(ni1)-gx(ni1+1)
gz(nk1-n)=gz(nk1-n+1)+gz(nk1)-gz(nk1+1)
gx(NTPX-2+n)=gx(NTPX-2+n-1)+gx(NTPX-2)-gx(NTPX-2-1)
gz(NTPZ-2+n)=gz(NTPZ-2+n-1)+gz(NTPZ-2)-gz(NTPZ-2-1)
end do

do k=1,NTPZ-1
gzhalf(k)=gz(k)+(gz(k+1)-gz(k))/2.0
end do
gzhalf(NTPZ)=gz(NTPZ)+(gz(NTPZ)-gz(NTPZ-1))/2.0

gtopoz=gzhalf(NTPZ-2)
topogrid=gzhalf(NTPZ-2)

!-----------------------------------------------------------------------------
print *, "output grids coordinate ..."

!output grid
filenm=grid_coordfnm_get()
call coord_create(filenm)

!-----------------------------------------------------------------------------

call grid_destroy
call dealloc_local

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

!*************************************************************************
!*                    PART-I  alloc and dealloc                          *
!*************************************************************************
subroutine alloc_local
  integer :: ierr
  allocate(indxVx(SEIS_GEO,ni*nk),stat=ierr); indxVx=0
  allocate(indxVz(SEIS_GEO,ni*nk),stat=ierr); indxVz=0
  allocate(NzVx(  ni*nk),stat=ierr); NzVx=0.0
  allocate(NxVz(2,ni*nk),stat=ierr); NxVz=0.0
  allocate(gx(NTPX)); gx=0.0
  allocate(gz(NTPZ)); gz=0.0
  allocate(gzhalf(NTPZ)); gzhalf=0.0
  allocate(topogrid(nx)); topogrid=0.0
  allocate(gtopoz(NTPX)); gtopoz=0.0
  allocate(gfindx(NTPX)); gfindx=0
  allocate(gbit(NTPX)); gbit=.false.
  allocate(gbit_mask(NTPX)); gbit_mask=.false.
end subroutine alloc_local

subroutine dealloc_local
  if (allocated(gfindx)) deallocate(gfindx)
end subroutine dealloc_local

subroutine init_grid(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  character (len=SEIS_STRLEN) :: filenm
  integer :: fid,gid,i,k,ierr

  fid=1001
  allocate(G%x(NTPI)); G%x=0.0_SP
  allocate(G%z(NTPK)); G%z=0.0_SP

  open(fid,file=trim(fnm_conf),status="old")
  gid=1002

  ! read in x coordinate
  call string_conf(fid,1,'coordx_filename',2,filenm)
  open(gid,file=trim(filenm),status='old')
  read(gid,*,iostat=ierr) i, stepx
  if (ierr<0) then
    write(*,*) trim(filenm)
    write(*,*) i,stepx
    call error_except('coordx: first line should be total number and stepx')
  end if
  if (i<NTPI) call error_except('no enough x grid points')
  if (i>NTPI) call warning_print('x grid in file is more than requaired')
  read(gid,*) ( G%x(i),i=1,NTPI )
  close(gid)

  ! read in z coordinate
  call string_conf(fid,1,'coordz_filename',2,filenm)
  open(gid,file=trim(filenm),status='old')
  read(gid,*,iostat=ierr) k, stepz
  if (ierr<0) then
    write(*,*) trim(filenm)
    write(*,*) k,stepz
    call error_except('coordz: first line should be total number and stepz')
  end if
  if (k<NTPK) call error_except('no enough z grid points')
  if (k>NTPK) call warning_print('z grid in file is more than requaired')
  read(gid,*) ( G%z(k),k=1,NTPK )
  close(gid)

  close(fid)
end subroutine init_grid

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

!----------------------------------------------------------------
subroutine coord_create(filenm)
  character (len=*),intent(in) :: filenm

  integer ncid,ierr,oldMode
  integer xid,zid
  integer coordxid,coordzid,zhalfid,zgridid,topoid

  ierr=nf90_create( path=trim(filenm),cmode=nf90_clobber,ncid=ncid )
  ierr=nf90_set_fill(ncid,nf90_nofill,oldMode)
  ! -- define dim
  ierr=nf90_def_dim(ncid, 'I', nx, xid);     call nfseis_handle_err(ierr)
  ierr=nf90_def_dim(ncid, 'K', nz, zid);     call nfseis_handle_err(ierr)
  ! -- variable
  ierr=nf90_def_var(ncid, 'x', nf90_float, (/ xid /), coordxid )
  ierr=nf90_def_var(ncid, 'z', nf90_float, (/ zid /), coordzid )
  ierr=nf90_def_var(ncid, 'z_of_force', nf90_float, (/ zid /), zhalfid )
  ierr=nf90_def_var(ncid, 'topo_on_grid', nf90_float, (/ xid /), zgridid )
  ierr=nf90_def_var(ncid, 'topography', nf90_float, (/ xid /), topoid )
  ierr=nf90_put_att(ncid,NF90_GLOBAL,'stepx',stepx)
  ierr=nf90_put_att(ncid,NF90_GLOBAL,'stepz',stepz)
  !--
  ierr=nf90_enddef(ncid)
  !--
  ierr=nf90_put_var(ncid,coordxid,gx(ngx1:ngx2),(/1/),(/nx/),(/1/))
  ierr=nf90_put_var(ncid,coordzid,gz(ngz1:ngz2),(/1/),(/nz/),(/1/))
  ierr=nf90_put_var(ncid,zhalfid,gzhalf(ngz1:ngz2),(/1/),(/nz/),(/1/))
  ierr=nf90_put_var(ncid,zgridid,topogrid,(/1/),(/nx/),(/1/))
  ierr=nf90_put_var(ncid,topoid,gtopoz(ngx1:ngx2),(/1/),(/nx/),(/1/))
  ierr=nf90_close(ncid); call nfseis_handle_err(ierr)
end subroutine coord_create

subroutine error_except(msg)
  character (len=*),intent(in) :: msg
  print *, trim(msg)
  stop 1
end subroutine error_except
subroutine warning_print(msg)
  character (len=*),intent(in) :: msg
  print *, " warning :"//trim(msg)
end subroutine warning_print

end program seis2d_grid

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
