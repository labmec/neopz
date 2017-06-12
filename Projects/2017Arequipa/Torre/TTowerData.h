#ifndef TTowerDataH
#define TTowerDataH
#include "json.h"
#include <string>
#include <map>

/** Tower structural data
 * @author Tiago Forti, IST Institute, CEATI project
 * @since March 13, 2017
 */
class TTowerData
{

    public:
        void Write(std::string path);
        void Read(std::string path);

    public:

        struct TMatProp
        {
                double fRho, fE, fG;
        };

        struct TSecProp
        {
                double fA, fJt, fIy, fIz;
                std::string fMaterial;
        };

        struct TJoint
        {
                double x, y, z;
        };

        struct TJointRestaint
        {
                int fU1, fU2, fU3, fR1, fR2, fR3;
        };

        struct TIncid
        {
                int fJointI, fJointJ;
        };

        /** Material properties: MaterialId -> data
         * TABLE:  "MATERIAL PROPERTIES 02 - BASIC MECHANICAL PROPERTIES"
         * Material=A992Fy50   UnitWeight=76972.8639422648   UnitMass=7849.04737995992   E1=199947978795.958   G12=76903068767.676
         *          Keyword                                  fRho                        fE                    fG
         */
        std::map<std::string, TMatProp> fMaterialProp;

        /** Section properties: SectionId -> data
         * TABLE:  "FRAME SECTION PROPERTIES 01 - GENERAL"
         * SectionName=FSEC1   Material=A992Fy50   Area=0.0055040276   TorsConst=1.13171949538729E-07   I33=2.06162439662021E-04 I22=3.30542221844342E-06
         *             Keyword fMaterial           fA                  fJt                              fIy                      fIz
         */
        std::map<std::string, TSecProp> fSectionProp;

        /** Joints
         * TABLE:  "JOINT COORDINATES"
         * Joint=1   CoordSys=GLOBAL   CoordType=Cartesian   XorR=0   Y=0   Z=0   SpecialJt=No   GlobalX=0   GlobalY=0   GlobalZ=0
         * Keyword                                           x        y     z
         */
        std::map<int, TJoint> fJoints;

        /* Frames
         * TABLE:  "CONNECTIVITY - FRAME"
         * Frame=1   JointI=1   JointJ=2   IsCurved=No   Length=17.3205080756888   CentroidX=5   CentroidY=5   CentroidZ=5
         * Keyword   fJointI    fJointJ
         */
        std::map<int, TIncid> fFrameIncid;

        /* Frame assignment - section
         * TABLE:  "FRAME SECTION ASSIGNMENTS"
         * Frame=2   SectionType=Angle   AutoSelect=N.A.   AnalSect=L90X90X7   DesignSect=L90X90X7   MatProp=Default
         * Keyword                                         string
         */
        std::map<int, std::string> fFrameSection;

        /* Frame assignement - section rotation
         * TABLE:  "FRAME LOCAL AXES ASSIGNMENTS 1 - TYPICAL"
         * Frame=38   Angle=180   AdvanceAxes=No
         * Keyword    double
         */
        std::map<int, double> fFrameSectionRotation;

        /** Supports
         * TABLE:  "JOINT RESTRAINT ASSIGNMENTS"
         * Joint=1   U1=Yes   U2=Yes   U3=Yes   R1=Yes   R2=Yes   R3=Yes
         * Keyword   fU1      fU2      fU3      fR1      fR2      fR3
         */
        std::map<int, TJointRestaint> fJointRestr;

};

class JsonParser
{
    public:
        JsonParser(const std::string inputPath);
        ~JsonParser();
        void Read(TTowerData & ttData);
        void Write(TTowerData & ttData);
    private:
        void DefineTables();
        // Read operations
        void ReadMaterialProp(TTowerData & ttData);
        void ReadSectionProp(TTowerData & ttData);
        void ReadJoints(TTowerData & ttData);
        void ReadFrameIncid(TTowerData & ttData);
        void ReadFrameSection(TTowerData & ttData);
        void ReadFrameSectionRotation(TTowerData & ttData);
        void ReadJointRestr(TTowerData & ttData);
        // Write operations
        void WriteMaterialProp(TTowerData & ttData);
        void WriteSectionProp(TTowerData & ttData);
        void WriteJoints(TTowerData & ttData);
        void WriteFrameIncid(TTowerData & ttData);
        void WriteFrameSection(TTowerData & ttData);
        void WriteFrameSectionRotation(TTowerData & ttData);
        void WriteJointRestr(TTowerData & ttData);
    private:
        const char * MATERIAL_PROP_TABLE;
        const char * SECTION_PROP_TABLE;
        const char * JOINTS_TABLE;
        const char * FRAME_INCID_TABLE;
        const char * FRAME_SECTION_TABLE;
        const char * FRAME_SECTION_ROTATION_TABLE;
        const char * JOINT_RESTR_TABLE;
        std::string fPath;
        Json::Value fRoot;

};
#endif

