//
//  TestGenGrid.cpp
//  PZ
//
//  Created by Pedro Lima on 21/jun/2021.
//

#include "pzgmesh.h"
#include "TPZGenGrid3D.h"
#include "TPZGenGrid2D.h"

#include <stdio.h>
#include <iostream>

#include <catch2/catch.hpp>

// Forward declaring
namespace gengridtest{
    /** @brief Testing aspects of TPZGenGrid3D*/
    void TestGenGrid3D(MMeshType meshtype);
    /** @brief Testing aspects of TPZGenGrid2D*/
    void TestGenGrid2D(MMeshType meshtype);
    /** @brief Require positive detjacobian of all elements in the mesh*/
    void TestJacobian(TPZGeoMesh* gmesh);
}

    /** @brief Test case to verify consistency of TPZGenGrid3D*/
    TEST_CASE("GenGrid3D_tests_1","[GenGrid3D_tests]")
    {
        REQUIRE_THROWS(gengridtest::TestGenGrid3D(MMeshType::EQuadrilateral));  /*For 2D meshes consisting of quadrilaterals*/
        REQUIRE_THROWS(gengridtest::TestGenGrid3D(MMeshType::ETriangular));     /*For 2D meshes where each quadrilateral is split into two triangles*/
        gengridtest::TestGenGrid3D(MMeshType::EHexahedral);     /*For 3D meshes consisting of hexahedra*/
        gengridtest::TestGenGrid3D(MMeshType::ETetrahedral);    /*For 3D meshes where each hexahedron is divided into five tetrahedra*/
        gengridtest::TestGenGrid3D(MMeshType::EPyramidal);      /*For 3D meshes where each hexahedron is divided into six pyramids*/
        gengridtest::TestGenGrid3D(MMeshType::EPrismatic);      /*For 3D meshes where each hexahedron is divided into two prisms*/
        gengridtest::TestGenGrid3D(MMeshType::EHexaPyrMixed);   /*For 3D meshes where alternating cubes are divided into six pyramids*/
        REQUIRE_THROWS(gengridtest::TestGenGrid3D(MMeshType::ENoType));             
    }
/** @brief Test case to verify consistency of TPZGenGrid2D*/
    TEST_CASE("GenGrid2D_tests_1","[GenGrid2D_tests]")
    {
        gengridtest::TestGenGrid2D(MMeshType::EQuadrilateral);  /*For 2D meshes consisting of quadrilaterals*/
        gengridtest::TestGenGrid2D(MMeshType::ETriangular);     /*For 2D meshes where each quadrilateral is split into two triangles*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::EHexahedral));     /*For 2D meshes consisting of hexahedra*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::ETetrahedral));    /*For 2D meshes where each hexahedron is divided into five tetrahedra*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::EPyramidal));      /*For 2D meshes where each hexahedron is divided into six pyramids*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::EPrismatic));      /*For 2D meshes where each hexahedron is divided into two prisms*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::EHexaPyrMixed));   /*For 2D meshes where alternating cubes are divided into six pyramids*/
        REQUIRE_THROWS(gengridtest::TestGenGrid2D(MMeshType::ENoType));             
    }
    /** Suggestions for future tests to add:
     * - Initialize and test ALL grid generating tools available in PZ
    */


namespace gengridtest{


    /** @brief Testing aspects of TPZGenGrid3D*/
    void TestGenGrid3D(MMeshType meshtype){
        // Create GenGrid
        TPZGenGrid3D gengrid({-6,-6,-6},{6,6,6},{3,3,3},meshtype);
        TPZGeoMesh* gmesh = gengrid.BuildVolumetricElements(1);
        // Check determinant of jacobian matrix for all elements in the mesh
        TestJacobian(gmesh);

        delete gmesh;
    }

    /** @brief Testing aspects of TPZGenGrid3D*/
    void TestGenGrid2D(MMeshType meshtype){
        // Create GenGrid
        TPZGenGrid2D gengrid({3,3},{-2.,-2.},{2.,2.});
        gengrid.SetElementType(meshtype);
        TPZGeoMesh* gmesh = new TPZGeoMesh;
        gengrid.Read(gmesh,1);
        // Check determinant of jacobian matrix for all elements in the mesh
        TestJacobian(gmesh);

        delete gmesh;
    }







    /** @brief Require positive detjacobian of all elements in the mesh*/
    void TestJacobian(TPZGeoMesh* gmesh){

        constexpr REAL notcomputed = -1024.0;
        static std::string testName = __PRETTY_FUNCTION__;

        for(auto gel : gmesh->ElementVec()){
            if(!gel) continue;

            // Compute jacobian of geomtric map
            TPZVec<REAL> qsi(gel->Dimension(),0.); TPZFMatrix<double> jac; TPZFMatrix<double> axes; REAL detjac = notcomputed; TPZFMatrix<double> jacinv;
            gel->Jacobian(qsi,jac,axes,detjac,jacinv);

            // Require that det(J) be positive
            const REAL zero = ZeroTolerance();
            const bool cond = detjac > zero;
            CAPTURE(gel->Index(),jac,detjac);
            REQUIRE(cond);
        }
    }
}//namespace gengridtest
