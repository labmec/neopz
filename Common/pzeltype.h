/**
 * @file
 * @brief Defines enum MElementType and contains the implementation of MElementType_NNodes(...) functions
 */
#ifndef PZELTYPEH
#define PZELTYPEH

#include <string>
#include <iostream>
#include <cstdlib>
#include "pzerror.h"

/**
 * \addtogroup common
 * @{
 */

/**
 * Zero-dimensional: EPoint. \n
 * One-dimensional: EOned (element) EInterfacePoint (interface). \n
 * Two-dimensional: ETriangle EQuadrilateral (element) EInterfaceLinear (interface). \n
 * Three-dimensional: ETetrahedro EPiramide EPrisma ECube (element) EInterfaceSurface (interface). \n
 * n-dimensional: ESubstructure EGlobLoc EDiscontinuous EInterfaceDisc EAgglomerate
 *
 * @brief Define the element types.
 */

// $Id: pzeltype.h,v 1.9 2011-04-05 19:32:54 calle Exp $
/**
   @enum MElementType
 * Defines the element types
 * @param EPoint            element 0D - type point        -  associated index 0
 * @param EOned             element 1D - type oned         -  associated index 1
 * @param ETriangle         element 2D - type triangle     -  associated index 2
 * @param EQuadrilateral    element 2D - type quad         -  associated index 3
 * @param ETetraedro        element 3D - type tetraedro    -  associated index 4
 * @param EPiramide         element 3D - type piramide     -  associated index 5
 * @param EPrisma           element 3D - type prisma       -  associated index 6
 * @param ECube             element 3D - type cube         -  associated index 7
 * @param EPolygonal        element ?? - type ??           -  associated index 8
 * @param EInterface        element nD - type interface    -  associated index 9
 * @param EInterfacePoint   element 1D - type interface    -  associated index 10
 * @param EInterfaceLinear  element 2D - type interface    -  associated index 11
 * @param EInterfaceSurface element 3D - type interface    -  associated index 12
 * @param ESubstructure     element nD - type submesh      -  associated index 13
 * @param EGlobLoc          element nD - type global-local -  associated index 14
 * @param EDiscontinuous    element nD - type discontinuous-  associated index 15
 * @param EInterfaceDisc    element nD - type interface    -  associated index 16
 * @param EAgglomerate      element nD - type agglomerate  -  associated index 17
 * @param ENoType           element 0D - type none         -  associated index 18
 */

enum MElementType
{
	/*0*/	EPoint,
	/*1*/	EOned,
	/*2*/	ETriangle,
	/*3*/	EQuadrilateral,
	/*4*/	ETetraedro,
	/*5*/	EPiramide,
	/*6*/	EPrisma,
	/*7*/	ECube,
	/*8*/	EPolygonal,
	/*9*/	EInterface,
	/*10*/	EInterfacePoint,
	/*11*/	EInterfaceLinear,
	/*12*/	EInterfaceSurface,
	/*13*/	ESubstructure,
	/*14*/	EGlobLoc,
	/*15*/	EDiscontinuous,
	/*16*/	EAgglomerate,
	/*17*/	ENoType,
	/*18*/	EInterfaceDisc = EInterface
};

/**
 * @brief Returns the number of nodes according to the type of the element
 */
inline int MElementType_NNodes(MElementType elType)
{
	switch(elType)
	{
        case(0)://point
        {
            return 1;
        }
		case(1)://line
		{
			return 2;
		}
		case(2)://triangle
		{
			return 3;
		}
		case(3)://quadrilateral
		{
			return 4;
		}
		case(4)://tetraedron
		{
			return 4;
		}
		case(5)://pyramid
		{
			return 5;
		}
		case(6)://prism
		{
			return 6;
		}
		case(7)://cube
		{
			return 8;
		}
		default:
		{
			std::cout << "ElementType not found!";
			DebugStop();
		}
	}
	DebugStop();
	return -1;
}

/**
 * @brief Returns the name of the element type.
 */
inline std::string MElementType_Name(MElementType elType)
{
	int elTypeId = elType;
	switch (elTypeId)
	{
		case 0:
		{
			return "EPoint";
		}
		case 1:
		{
			return "EOned";
		}
		case 2:
		{
			return "ETriangle";
		}
		case 3:
		{
			return "EQuadrilateral";
		}
		case 4:
		{
			return "ETetraedro";
		}
		case 5:
		{
			return "EPiramide";
		}
		case 6:
		{
			return "EPrisma";
		}
		case 7:
		{
			return "ECube";
		}
		case 8:
		{
			return "EPolygonal";
		}
		case 9:
		{
			return "EInterface";
		}
		case 10:
		{
			return "EInterfacePoint";
		}
		case 11:
		{
			return "EInterfaceLinear";
		}
		case 12:
		{
			return "EInterfaceSurface";
		}
		case 13:
		{
			return "ESubstructure";
		}
		case 14:
		{
			return "EGlobLoc";
		}
		case 15:
		{
			return "EDiscontinuous";
		}
		case 16:
		{
			return "EAgglomerate";
		}
		case 17:
		{
			return "ENoType";
		}
		case 18:
		{
			return "EInterfaceDisc = EInterface";
		}
		default:
		{
			return "ElementType not found!";
		}
	}
	DebugStop();
	return "";
}

/**
 * @}
 */

#endif
