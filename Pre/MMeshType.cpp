
#include <iostream>
#include "MMeshType.h"


std::ostream &operator<<(std::ostream &out, const MMeshType meshType)
{
    out << MMeshType_Name(meshType);
    return out;
}

std::string MMeshType_Name(const MMeshType meshType){
    switch(meshType){
        case MMeshType::ETriangular:
            return "Triangular";
        case MMeshType::EQuadrilateral:
            return "Quadrilateral";
        case MMeshType::ETetrahedral:
            return "Tetrahedral";
        case MMeshType::EPyramidal:
            return "Pyramidal";
        case MMeshType::EPrismatic:
            return "Prismatic";
        case MMeshType::EHexahedral:
            return "Hexahedral";
        case MMeshType::EHexaPyrMixed:
            return "Hex/Pyr";
        case MMeshType::ENoType:
            return "No Type";
    }
    DebugStop();
    return "No Type";
}

int MMeshType_Dimension(const MMeshType meshType){
    switch(meshType){
        case MMeshType::ETriangular:
        case MMeshType::EQuadrilateral:
            return 2;
        case MMeshType::ETetrahedral:
        case MMeshType::EPyramidal:
        case MMeshType::EPrismatic:
        case MMeshType::EHexahedral:
        case MMeshType::EHexaPyrMixed:
            return 3;
        case MMeshType::ENoType:
            return -1;
    }
    DebugStop();
    return -1;

}