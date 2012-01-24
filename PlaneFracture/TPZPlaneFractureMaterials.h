//
//  TPZPlaneFractureMaterials.h
//  PZ
//
//  Created by Cesar Lucci on 21/11/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//


/** @brief Material index of auxiliar 1D elements in 2D refpatterns generation */
const int __aux1DEl_Mat = -1;

//BCs plane fracture
const int __1DcrackTipMat = -10;//1D crack tip elements

const int __1DbulletMat = -20;//1D elements that introduces injection flow rate

const int __2DfractureMat_outside = -30;//2D elements outside fracture (i.e.: not fractured yet)

const int __2DfractureMat_inside = -40;//2D elements inside fracture (i.e.: already fractured)

//BCs surrounding
const int __2DfarfieldXZMat = -50;//2D elements in plane x,z (y_max > 0 = farfield)

const int __2Dleft_rightMat = -60;//2D elements in plane y,z (one in x=0 "left farfield" and another in x_max > 0 "right farfield")

const int __2Dtop_bottomMat = -70;//2D elements in plane x,y (vertical confinement)

//Domain
const int __3DrockMat_linear = 10;//3D elements that surround fracture plane (i.e.: porous media)

const int __3DrockMat_quarterPoint = 20;//3D elements surrounding crack tip