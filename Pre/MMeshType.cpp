
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



MMeshType StringToMeshtype(const std::string& name){
	if(name[0] != 'E'){PZError << "\nYou forgot the E";}
    else if(name == "EQuadrilateral"){	return MMeshType::EQuadrilateral;}
    else if(name == "ETriangular"){		return MMeshType::ETriangular;}
    else if(name == "EHexahedral"){		return MMeshType::EHexahedral;}
    else if(name == "ETetrahedral"){	return MMeshType::ETetrahedral;}
    else if(name == "EPyramidal"){		return MMeshType::EPyramidal;}
    else if(name == "EPrismatic"){		return MMeshType::EPrismatic;}
    else if(name == "EHexaPyrMixed"){	return MMeshType::EHexaPyrMixed;}
    else if(name == "ENoType"){			return MMeshType::ENoType;}
	
    PZError << "\nTried to interpret unrecognized MMeshType : \'" << name << "\'\n";
    DebugStop();
	return MMeshType::ENoType; // dummy return to remove warning
}