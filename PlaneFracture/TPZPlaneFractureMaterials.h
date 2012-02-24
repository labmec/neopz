//
//  TPZPlaneFractureMaterials.h
//  PZ
//
//  Created by Cesar Lucci on 21/11/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


/** @brief Material index of auxiliar 1D elements in 2D refpatterns generation */
const int __aux1DEl_Mat = -1;
const int __aux0DEl_Mat1 = -2;
const int __aux0DEl_Mat2 = -3;
const int __aux0DEl_Mat3 = -4;
const int __aux0DEl_Mat4 = -5;

//BCs plane fracture
const int __1DcrackTipMat = -10;//1D crack tip elements

const int __1DbulletMat = -20;//1D elements that introduces injection flow rate

const int __2DfractureMat_outside = -30;//2D elements outside fracture (i.e.: not fractured yet)

const int __2DfractureMat_inside = -40;//2D elements inside fracture (i.e.: already fractured)

//BCs surrounding
const int __2DfarfieldMat = -50;//2D elements in plane x,z (y_max > 0 = farfield)

const int __2DleftMat = -60;//2D elements in plane y,z (in x=0)

const int __2DrightMat = -70;//2D elements in plane y,z (in x=xmax)

const int __2DtopMat = -80;//2D elements in plane x,y (top confinement)

const int __2DbottomMat = -90;//2D elements in plane x,y (bottom confinement)

//Domain
const int __3DrockMat_linear = 10;//3D elements that surround fracture plane (i.e.: porous media)

const int __3DrockMat_quarterPoint = 20;//3D elements surrounding crack tip