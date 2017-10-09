/*---------------------------------------------------------------------------
! This header file contains the 4th-order Finite-Difference stencil.
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!--------------------------------------------------------------------------*/

#define DEFFDWET real,parameter,private :: MCC0=1.125, MCC1=0.0416666667
#define DEFFDWETPRO real,parameter :: MCC0=1.125, MCC1=0.0416666667
#define LenFD 2
#define LenFDS 1
#define LenFDL 2

/* #ifdef DataTypeDouble
#define SEISNC_DATATYPE NF90_DOUBLE
#define SEISMPI_DATATYPE MPI_DOUBLE_PRECISION
#else
#define SEISNC_DATATYPE NF90_FLOAT
#define SEISMPI_DATATYPE MPI_REAL
#endif */

/*********************************************************
 * ! fd for 3d     *
 *********************************************************/
#define m2d_FDxF(var,i,k)  ( MCC0*(var(i+1,k)-var(i,k)) \
-MCC1*(var(i+2,k)-var(i-1,k)) )/stepx
#define m2d_FDxB(var,i,k)  ( MCC0*(var(i,k)-var(i-1,k)) \
-MCC1*(var(i+1,k)-var(i-2,k)) )/stepx

#define m2d_FDzF(var,i,k)  ( MCC0*(var(i,k+1)-var(i,k)) \
-MCC1*(var(i,k+2)-var(i,k-1)) )/stepz
#define m2d_FDzB(var,i,k)  ( MCC0*(var(i,k)-var(i,k-1)) \
-MCC1*(var(i,k+1)-var(i,k-2)) )/stepz

/*********************************************************
 * ! fd for AFDA     *
 *********************************************************/
#define AFDA_HOCF_1DH(var,dev,i,k)  ( -stepz/22.*dev +577./528.*var(i,k+1) \
-201./176.*var(i,k)+9./176.*var(i,k-1)-1./528.*var(i,k-2) )/stepz
#define AFDA_HOCB_1DH(var,dev,i,k)  ( -stepz/22.*dev +577./528.*var(i,k) \
-201./176.*var(i,k-1)+9./176.*var(i,k-2)-1./528.*var(i,k-3) )/stepz

#define AFDA_DzTF_00DH(var,i,k)  ( -35./8.*var(i,k) \
+35./24.*var(i,k-1)-21./40.*var(i,k-2)+5./56.*var(i,k-3) )/stepz
#define AFDA_DzTB_00DH(var,i,k)  ( -35./8.*var(i,k-1) \
+35./24.*var(i,k-2)-21./40.*var(i,k-3)+5./56.*var(i,k-4) )/stepz

#define AFDA_DzTF_10DH(var,i,k)  ( 31./24.*var(i,k+1) \
-29./24.*var(i,k)+3./40.*var(i,k-1)-1./168.*var(i,k-2) )/stepz
#define AFDA_DzTB_10DH(var,i,k)  ( 31./24.*var(i,k) \
-29./24.*var(i,k-1)+3./40.*var(i,k-2)-1./168.*var(i,k-3) )/stepz

#define AFDA_DzF_05DH(var,i,k)  ( 11./12.*var(i,k+1) \
-17./24.*var(i,k)-3./8.*var(i,k-1)+5./24.*var(i,k-2)-1./24.*var(i,k-3) )/stepz
#define AFDA_DzB_05DH(var,i,k)  ( 11./12.*var(i,k) \
-17./24.*var(i,k-1)-3./8.*var(i,k-2)+5./24.*var(i,k-3)-1./24.*var(i,k-4) )/stepz

/* 
#define AFDA_DzF_05DH(var,i,k)  ( 23./24.*var(i,k+1) \
-7./8.*var(i,k)-1./8.*var(i,k-1)+1./24.*var(i,k-2) )/stepz
#define AFDA_DzB_05DH(var,i,k)  ( 23./24.*var(i,k) \
-7./8.*var(i,k-1)-1./8.*var(i,k-2)+1./24.*var(i,k-3) )/stepz
*/

/**************************************
 * ! fd for 1d *
 **************************************/
#define vec_FDxF(var,i)  ( MCC0*(var(i+1)-var(i)) \
-MCC1*(var(i+2)-var(i-1)) )/stepx
#define vec_FDxB(var,i)  ( MCC0*(var(i)-var(i-1)) \
-MCC1*(var(i+1)-var(i-2)) )/stepx

#define vec_FDzF(var,i)  ( MCC0*(var(i+1)-var(i)) \
-MCC1*(var(i+2)-var(i-1)) )/stepz
#define vec_FDzB(var,i)  ( MCC0*(var(i)-var(i-1)) \
-MCC1*(var(i+1)-var(i-2)) )/stepz

