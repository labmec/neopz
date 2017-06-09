#include "TTowerData.h"
//#include "TSWXDebug.h"
#include <fstream>
#include <iostream>

JsonParser::JsonParser(const std::string inputPath)
{
    fPath = inputPath;
    DefineTables();
}

JsonParser::~JsonParser()
{
    fRoot.clear();
}

void JsonParser::DefineTables()
{
    MATERIAL_PROP_TABLE = "MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES";
    SECTION_PROP_TABLE = "FRAME SECTION PROPERTIES 01 - GENERAL";
    JOINTS_TABLE = "JOINT COORDINATES";
    FRAME_INCID_TABLE = "CONNECTIVITY - FRAME";
    FRAME_SECTION_TABLE = "FRAME SECTION ASSIGNMENTS";
    FRAME_SECTION_ROTATION_TABLE = "FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL";
    JOINT_RESTR_TABLE = "JOINT RESTRAINT ASSIGNMENTS";
}

void JsonParser::Read(TTowerData & ttData)
{
    fRoot.clear();

    std::ifstream inputFile(fPath.c_str(), std::ifstream::binary);
    if (inputFile.fail())
    {
        std::cout << "Couldn't open input file json" << std::endl;
        return;
    }
    inputFile >> fRoot;

    ReadMaterialProp(ttData);
    ReadSectionProp(ttData);
    ReadJoints(ttData);
    ReadFrameIncid(ttData);
    ReadFrameSection(ttData);
    ReadFrameSectionRotation(ttData);
    ReadJointRestr(ttData);
}

void JsonParser::Write(TTowerData & ttData)
{
    fRoot.clear();

    WriteMaterialProp(ttData);
    WriteSectionProp(ttData);
    WriteJoints(ttData);
    WriteFrameIncid(ttData);
    WriteFrameSection(ttData);
    WriteFrameSectionRotation(ttData);
    WriteJointRestr(ttData);

    Json::StyledStreamWriter writer("    ");
    std::ofstream output;
    output.open(fPath.c_str());
    writer.write(output, fRoot);
    output.close();
}

void JsonParser::ReadMaterialProp(TTowerData & ttData)
{
    Json::Value table = fRoot[MATERIAL_PROP_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        std::string keyword = tableData[i]["Material"].asString();

        TTowerData::TMatProp matProp;
        matProp.fRho = tableData[i]["UnitMass"].asDouble();
        matProp.fE = tableData[i]["E1"].asDouble();
        matProp.fG = tableData[i]["G12"].asDouble();

        ttData.fMaterialProp[keyword] = matProp;
    }
}

void JsonParser::WriteMaterialProp(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<std::string, TTowerData::TMatProp>::iterator it;
    for (it = ttData.fMaterialProp.begin(); it != ttData.fMaterialProp.end(); it++)
    {
        Json::Value item;
        item["Material"] = it->first;
        item["UnitMass"] = it->second.fRho;
        item["E1"] = it->second.fE;
        item["G12"] = it->second.fG;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[MATERIAL_PROP_TABLE] = table;
}

void JsonParser::ReadSectionProp(TTowerData & ttData)
{
    Json::Value table = fRoot[SECTION_PROP_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        std::string keyword = tableData[i]["SectionName"].asString();

        TTowerData::TSecProp secProp;
        secProp.fMaterial = tableData[i]["Material"].asString();
        secProp.fA = tableData[i]["Area"].asDouble();
        secProp.fJt = tableData[i]["TorsConst"].asDouble();
        secProp.fIy = tableData[i]["I33"].asDouble();
        secProp.fIz = tableData[i]["I22"].asDouble();

        ttData.fSectionProp[keyword] = secProp;
    }
}

void JsonParser::WriteSectionProp(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<std::string, TTowerData::TSecProp>::iterator it;
    for (it = ttData.fSectionProp.begin(); it != ttData.fSectionProp.end(); it++)
    {
        Json::Value item;
        item["SectionName"] = it->first;
        item["Material"] = it->second.fMaterial;
        item["Area"] = it->second.fA;
        item["TorsConst"] = it->second.fJt;
        item["I33"] = it->second.fIy;
        item["I22"] = it->second.fIz;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[SECTION_PROP_TABLE] = table;
}

void JsonParser::ReadJoints(TTowerData & ttData)
{
    Json::Value table = fRoot[JOINTS_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        int keyword = tableData[i]["Joint"].asInt();

        TTowerData::TJoint joint;
        joint.x = tableData[i]["XorR"].asDouble();
        joint.y = tableData[i]["Y"].asDouble();
        joint.z = tableData[i]["Z"].asDouble();

        ttData.fJoints[keyword] = joint;
    }
}

void JsonParser::WriteJoints(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<int, TTowerData::TJoint>::iterator it;
    for (it = ttData.fJoints.begin(); it != ttData.fJoints.end(); it++)
    {
        Json::Value item;
        item["Joint"] = it->first;
        item["XorR"] = it->second.x;
        item["Y"] = it->second.y;
        item["Z"] = it->second.z;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[JOINTS_TABLE] = table;
}

void JsonParser::ReadFrameIncid(TTowerData & ttData)
{
    Json::Value table = fRoot[FRAME_INCID_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        int keyword = tableData[i]["Frame"].asInt();

        TTowerData::TIncid frame;
        frame.fJointI = tableData[i]["JointI"].asInt();
        frame.fJointJ = tableData[i]["JointJ"].asInt();

        ttData.fFrameIncid[keyword] = frame;
    }
}

void JsonParser::WriteFrameIncid(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<int, TTowerData::TIncid>::iterator it;
    for (it = ttData.fFrameIncid.begin(); it != ttData.fFrameIncid.end(); it++)
    {
        Json::Value item;
        item["Frame"] = it->first;
        item["JointI"] = it->second.fJointI;
        item["JointJ"] = it->second.fJointJ;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[FRAME_INCID_TABLE] = table;
}

void JsonParser::ReadFrameSection(TTowerData & ttData)
{
    Json::Value table = fRoot[FRAME_SECTION_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        int keyword = tableData[i]["Frame"].asInt();

        ttData.fFrameSection[keyword] = tableData[i]["AnalSect"].asString();
    }
}

void JsonParser::WriteFrameSection(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<int, std::string>::iterator it;
    for (it = ttData.fFrameSection.begin(); it != ttData.fFrameSection.end(); it++)
    {
        Json::Value item;
        item["Frame"] = it->first;
        item["AnalSect"] = it->second;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[FRAME_SECTION_TABLE] = table;
}

void JsonParser::ReadFrameSectionRotation(TTowerData & ttData)
{
    if (fRoot.isMember(FRAME_SECTION_ROTATION_TABLE))
    {
        Json::Value table = fRoot[FRAME_SECTION_ROTATION_TABLE];
        Json::Value tableData = table["data"];

        int numData = tableData.size();
        for (int i = 0; i < numData; i++)
        {
            int keyword = tableData[i]["Frame"].asInt();

            ttData.fFrameSectionRotation[keyword] = tableData[i]["Angle"].asDouble();
        }
    }
}

void JsonParser::WriteFrameSectionRotation(TTowerData & ttData)
{
    int numData = ttData.fFrameSectionRotation.size();
    if (numData > 0)
    {
        Json::Value table;
        Json::Value tableData(Json::arrayValue);

        std::map<int, double>::iterator it;
        for (it = ttData.fFrameSectionRotation.begin(); it != ttData.fFrameSectionRotation.end(); it++)
        {
            Json::Value item;
            item["Frame"] = it->first;
            item["Angle"] = it->second;

            tableData.append(item);
        }
        table["data"] = tableData;

        fRoot[FRAME_SECTION_ROTATION_TABLE] = table;
    }
}

void JsonParser::ReadJointRestr(TTowerData & ttData)
{
    Json::Value table = fRoot[JOINT_RESTR_TABLE];
    Json::Value tableData = table["data"];

    int numData = tableData.size();
    for (int i = 0; i < numData; i++)
    {
        int keyword = tableData[i]["Joint"].asInt();

        TTowerData::TJointRestaint joint;
        joint.fU1 = tableData[i]["U1"].asInt();
        joint.fU2 = tableData[i]["U2"].asInt();
        joint.fU3 = tableData[i]["U3"].asInt();
        joint.fR1 = tableData[i]["R1"].asInt();
        joint.fR2 = tableData[i]["R2"].asInt();
        joint.fR3 = tableData[i]["R3"].asInt();

        ttData.fJointRestr[keyword] = joint;
    }
}

void JsonParser::WriteJointRestr(TTowerData & ttData)
{
    Json::Value table;
    Json::Value tableData(Json::arrayValue);

    std::map<int, TTowerData::TJointRestaint>::iterator it;
    for (it = ttData.fJointRestr.begin(); it != ttData.fJointRestr.end(); it++)
    {
        Json::Value item;
        item["Joint"] = it->first;
        item["U1"] = it->second.fU1;
        item["U2"] = it->second.fU2;
        item["U3"] = it->second.fU3;
        item["R1"] = it->second.fR1;
        item["R2"] = it->second.fR2;
        item["R3"] = it->second.fR3;

        tableData.append(item);
    }
    table["data"] = tableData;

    fRoot[JOINT_RESTR_TABLE] = table;
}

void TTowerData::Read(std::string path)
{
    JsonParser reader(path);
    reader.Read(*this);
}

void TTowerData::Write(std::string path)
{
    JsonParser writer(path);
    writer.Write(*this);
}
