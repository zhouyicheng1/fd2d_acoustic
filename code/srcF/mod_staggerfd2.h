/*---------------------------------------------------------------------------
! This header file contains the 2rd-order Finite-Difference stencil.
!
! Author: 
!       Wenzhong CAO    Email: caowz@mail.ustc.edu.cn
!       Wei ZHANG       Email: zhangwei.zw@gmail.com
! Copyright (C) 2017 Wenzhong CAO & Wei ZHANG
!--------------------------------------------------------------------------*/


#define DEFFDWET real,parameter,private :: MCC0=1.0
#define DEFFDWETPRO real,parameter :: MCC0=1.0
#define LenFD 1
#define LenFDS 1
#define LenFDL 1

/*********************************************************
 * ! fd for 3d     *
 *********************************************************/
#define m2d_FDxF(var,i,k)  ( MCC0*(var(i+1,k)-var(i,k)) )/stepx
#define m2d_FDxB(var,i,k)  ( MCC0*(var(i,k)-var(i-1,k)) )/stepx

#define m2d_FDzF(var,i,k)  ( MCC0*(var(i,k+1)-var(i,k)) )/stepz
#define m2d_FDzB(var,i,k)  ( MCC0*(var(i,k)-var(i,k-1)) )/stepz

/**************************************
 * ! fd for 1d *
 **************************************/
#define vec_FDxF(var,i)  ( MCC0*(var(i+1)-var(i)) )/stepx
#define vec_FDxB(var,i)  ( MCC0*(var(i)-var(i-1)) )/stepx

#define vec_FDzF(var,i)  ( MCC0*(var(i+1)-var(i)) )/stepz
#define vec_FDzB(var,i)  ( MCC0*(var(i)-var(i-1)) )/stepz

/*********************************************************
 * ! fd for AFDA     *
 *********************************************************/
/* this not used in 2nd-order, only for compile. */
#define AFDA_HOCF_1DH(var,dev,i,k)  0
#define AFDA_HOCB_1DH(var,dev,i,k)  0

#define AFDA_DzTF_00DH(var,i,k)     0
#define AFDA_DzTB_00DH(var,i,k)     0

#define AFDA_DzTF_10DH(var,i,k)     0 
#define AFDA_DzTB_10DH(var,i,k)     0

#define AFDA_DzF_05DH(var,i,k)      0 
#define AFDA_DzB_05DH(var,i,k)      0

