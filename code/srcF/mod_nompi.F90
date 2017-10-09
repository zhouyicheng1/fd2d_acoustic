module nompi_mod

!----------------------------------------------------------------------------
! This module is used for no mpi calculate
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
use para_mod
use string_mod

implicit none

private
public :: swmpi_init,         &
  cart_creat,                 &
  swmpi_reinit_para,          &
  swmpi_set_gindx,            &
  swmpi_globi,swmpi_globk,    &
  swmpi_locli,swmpi_loclk,    &
  point_in_thisnode,          &
  swmpi_time_init,            &
  swmpi_time_write,           &
  swmpi_time_end,             &
  swmpi_except

!-- #ifdef DataTypeDouble
!-- integer,parameter,public :: SEISMPI_DATATYPE=MPI_DOUBLE_PRECISION
!-- #else
!-- integer,parameter,public :: SEISMPI_DATATYPE=MPI_REAL
!-- #endif

integer,public :: ndims
integer,dimension(SEIS_GEO),public :: dims
integer,public :: SWMPI_COMM !communicator include topo information
integer,public :: myid
integer,dimension(SEIS_GEO),public :: thisid  !the thread coords in CART topo
integer,dimension(SEIS_GEO,2),public :: neigid

integer,public :: DTypeXS,DTypeXL,DTypeZS,DTypeZL

logical,public :: freenode,masternode
logical,dimension(SEIS_GEO,2),public :: absnode
character (len=SEIS_STRLEN),public :: fnm_swmpi

logical,dimension(SEIS_GEO) :: periods
logical :: reorder

character (len=8) :: str_d0,str_d1,str_d2
character (len=10):: str_t0,str_t1,str_t2
real (kind=8) :: wtime0,wtime1,wtime2
integer :: wtimen0

!-----------------------------------------------------------------------!
contains
!-----------------------------------------------------------------------!

subroutine swmpi_init(fnm_conf)
  character (len=*),intent(in) :: fnm_conf
  integer fid,n
  fid=1001
  open(fid,file=trim(fnm_conf),status="old")
  !call string_conf(fid,1,'ndims',2,ndims)
  ndims=SEIS_GEO
  do n=1,ndims
    call string_conf(fid,1,'dims',n+1,dims(n))
  end do
  periods=.false.
  reorder=.true.
  close(fid)
end subroutine swmpi_init
subroutine cart_creat
  integer m,n
  myid=0; thisid=0
  freenode = .true.
  absnode  = .true.

  masternode=.false.
  if (myid==0) masternode=.true.

end subroutine cart_creat
subroutine swmpi_reinit_para
  ! init ngi1 etc
  ngi1=swmpi_globi(ni1,thisid(1))
  ngi2=swmpi_globi(ni2,thisid(1))
  ngk1=swmpi_globk(nk1,thisid(2))
  ngk2=swmpi_globk(nk2,thisid(2))
  ngx1=swmpi_globi(nx1,thisid(1))
  ngx2=swmpi_globi(nx2,thisid(1))
  ngz1=swmpi_globk(nz1,thisid(2))
  ngz2=swmpi_globk(nz2,thisid(2))
  point_in_this=(/ ngi1,ngi2,ngk1,ngk2 /)
  if (thisid(1)==0)        point_in_this(1)=ngx1
  if (thisid(1)==dims(1)-1) point_in_this(2)=ngx2
  if (thisid(2)==0)        point_in_this(3)=ngz1
  if (thisid(2)==dims(2)-1) point_in_this(4)=ngz2

  NTPI=ni*dims(1); NTPK=nk*dims(2)
  NTPX=NTPI+(nx-ni); NTPZ=NTPK+(nz-nk)
  npi1=thisid(1)*ni+1; npi2=(thisid(1)+1)*ni
  npk1=thisid(2)*nk+1; npk2=(thisid(2)+1)*nk
end subroutine swmpi_reinit_para

subroutine swmpi_set_gindx(n_i,n_k)
  integer,intent(in) :: n_i,n_k
  ! init ngi1 etc
  ngi1=swmpi_globi(ni1,n_i)
  ngi2=swmpi_globi(ni2,n_i)
  ngk1=swmpi_globk(nk1,n_k)
  ngk2=swmpi_globk(nk2,n_k)
  ngx1=swmpi_globi(nx1,n_i)
  ngx2=swmpi_globi(nx2,n_i)
  ngz1=swmpi_globk(nz1,n_k)
  ngz2=swmpi_globk(nz2,n_k)
  point_in_this=(/ ngi1,ngi2,ngk1,ngk2 /)
  if (n_i==0)         point_in_this(1)=ngx1
  if (n_i==dims(1)-1) point_in_this(2)=ngx2
  if (n_k==0)         point_in_this(3)=ngz1
  if (n_k==dims(2)-1) point_in_this(4)=ngz2
  NTPI=ni*dims(1); NTPK=nk*dims(2)
  NTPX=NTPI+(nx-ni); NTPZ=NTPK+(nz-nk)
  npi1=n_i*ni+1; npi2=(n_i+1)*ni
  npk1=n_k*nk+1; npk2=(n_k+1)*nk
end subroutine swmpi_set_gindx

!*************************************************************************
!* nx,ny,nz, et al, should be inited before this subroutine              *
!*************************************************************************

subroutine swmpi_time_init(filenm,ntime)
  character (len=*),intent(in) :: filenm
  integer,intent(in) :: ntime
  integer fid
  fid=1009
  if (myid/=0) return
  call CPU_TIME(wtime0)
  call date_and_time(date=str_d0,time=str_t0)

  if (ntime==0) then
    open(fid,file=trim(filenm),status="unknown")
    write(fid,*) "# seis3d_wave run time log"
  else
    open(fid,file=trim(filenm),status="old",position="append")
    write(fid,*)
    write(fid,*) "# seis3d_wave restart from ntime=",ntime
  end if
  write(fid,*)
  write(fid,*) 'the program begins from ',str_d0,  &
    ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
  write(fid,*)
  !write(fid,*)  'each step uses time'
  !write(fid,'(4a10)') 'step','hour','minute','second'
  close(fid)
  wtime1=wtime0
  wtimen0=ntime
end subroutine swmpi_time_init
subroutine swmpi_time_write(ntime,filenm)
  character (len=*),intent(in) :: filenm
  integer,intent(in) :: ntime
  integer :: fid
  !integer h,m
  real (kind=8) :: s,s0,s1
  if (myid/=0) return
  fid=1009
  !-- wtime2=MPI_WTIME()
  call CPU_TIME(wtime2)
  s=wtime2-wtime1
  !h=int(s/3600)
  !m=int((s-h*3600)/60)
  !s=s-h*3600-m*60
  s0=(wtime2-wtime0)/3600.0
  s1=s0/(ntime-wtimen0)*(nt-wtimen0)
  open(fid,file=trim(filenm),status="old",position="append")
  !-- write(fid,'(i6,a10,g12.5,a4,g12.5,a17,g12.5,a11)') &
  !--         ntime,'step uses',s,'s,  ',                &
  !--         s0,' hours passed and', s1-s0, ' hours left'
  write(fid,'(i6,a15)') ntime,'step passed'
  close(fid)

  wtime1=wtime2
end subroutine swmpi_time_write
subroutine swmpi_time_end(filenm)
  character (len=*),intent(in) :: filenm
  integer fid
  integer d,h,m
  real (kind=8) :: s,wtimend

  if (myid/=0) return
  fid=1009
  !-- s=MPI_WTIME()-wtime0
  call CPU_TIME(wtimend)
  s=wtimend-wtime0
  d =int(s/3600.0/24.0)
  h =int((s-d*3600*24)/3600)
  m =int(s-d*3600*24-h*3600)/60
  call date_and_time(date=str_d1,time=str_t1)
  open(fid,file=trim(filenm),status="old",position="append")
  write(fid,*) '------------------------------------'
  write(fid,*) 'the program'
  write(fid,*) '  begins from ',str_d0,  &
    ',',str_t0(1:2),'/',str_t0(3:4),'/',str_t0(5:10)
  write(fid,*) '  finish at ',str_d1,  &
    ',',str_t1(1:2),'/',str_t1(3:4),'/',str_t1(5:10)
  !-- write(fid, * )  'the total run time is'
  !-- write(fid,'(3(i10,a10))') d,'day',h,'hour',m,'minute'
  close(fid)
end subroutine swmpi_time_end

function swmpi_globi(i,n_i) result(gi)
  integer,intent(in) :: i,n_i
  integer :: gi
  gi=(i-ni1+1)+n_i*ni+(ni1-1) ! term of (ni1-1) is for clarity
end function swmpi_globi
function swmpi_globk(k,n_k) result(gk)
  integer,intent(in) :: k,n_k
  integer :: gk
  gk=(k-nk1+1)+n_k*nk+(nk1-1)
end function swmpi_globk

function swmpi_locli(gi,n_i) result(i)
  integer,intent(in) :: gi,n_i
  integer :: i
  !i=mod(gi-(ni1-nx1),ni)+(ni1-nx1)
  i=gi-n_i*ni
end function swmpi_locli
function swmpi_loclk(gk,n_k) result(k)
  integer,intent(in) :: gk,n_k
  integer :: k
  !k=mod(gk-(nk1-nz1),nk)+(nk1-nz1)
  k=gk-n_k*nk
end function swmpi_loclk

function point_in_thisnode(i,k) result(iflag)
  integer,intent(in) :: i,k
  logical iflag
  integer glob_i1,glob_i2,glob_k1,glob_k2
  glob_i1=swmpi_globi(ni1,thisid(1)); glob_i2=swmpi_globi(ni2,thisid(1))
  glob_k1=swmpi_globk(nk1,thisid(1)); glob_k2=swmpi_globk(nk2,thisid(2))
  iflag=.false.
  if (      i>=glob_i1 .and. i<=glob_i2  &
    .and. k>=glob_k1 .and. k<=glob_k2 ) then
    iflag=.true.
  end if
end function point_in_thisnode

subroutine swmpi_except(msg)
  character (len=*),intent(in) :: msg
  integer :: ierr
  print *, trim(msg)
  stop 1
end subroutine swmpi_except

end module nompi_mod

! vim:ft=fortran:ts=2:sw=2:nu:et:ai:
