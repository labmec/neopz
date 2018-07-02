
#define _USE_MATH_DEFINES

#include <math.h>
#include <cmath>
#include "pzerror.h"

#include "TTowerData.h"

#include "pzelast3dGD.h"
#include "TPZEulerBernoulliBeamData.h"
#include "TPZEulerBernoulliBeam.h"
#include "TPZEulerBernoulliBC.h"
#include "TPZNewmarkAnalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

//#include "TSWXPowerEigenvalueMethod.h"
#include "TSWXExportFrameMesh.h"

/** Computes the first eigenmode of the structure.
 * It is implemented a modal analysis M.u" + K.u = 0
 * The eigenvalue problem is solved by the power algorithm
 * (see Mechanical Vibrations - Theory and Application to Structural Dynamics
 *      M. Géradin / D. Rixen - WILEY, 1994)
 *
 * @param Restraints indicates degrees of freedom with
 * displacement equal to ZERO.
 *
 * @author Tiago Forti
 * @author Prof. Philippe Devloo (phil@fec.unicamp.br)
 *
 * @since March 14, 2017
 */
class TSWXPowerEigenvalueMethod{

  public:

	  bool EigenModes(TPZAutoPointer<TPZMatrix<STATE> > Stiffness, TPZAutoPointer<TPZMatrix<STATE> > &Mass,
                  const TPZVec<int> & Restraints, REAL tol, int niter,
				  REAL & EigenValue, TPZFMatrix<STATE> & EigenVector);

};


TPZCompMesh * CreateCompMesh(TTowerData &data, TPZAutoPointer< TPZEulerBernoulliBeamData > & PropertyData,int porder);
void GetEquationsRestrained(TPZCompMesh & cmesh,TPZVec<int> &Restraints);

int main() {

  TPZMaterial::gBigNumber = 1e22;

  TPZAutoPointer< TPZEulerBernoulliBeamData > PropertyData;

  TTowerData data;
  std::string arq = PZSOURCEDIR;
//	std::string arq = "D:/SimulacaoNumerica/Simworx/VersionOnNeoPZ/NeoPZ";
	arq += "/Projects/2017Arequipa/Torre/torre.json";
  //data.Read(  "D:\\SimulacaoNumerica\\programaSimworx\\NeoPZ\\Projects\\Arequipa_Seminario\\Torre\\torre.json");
  data.Read(arq);
  int porder = 1;
  TPZCompMesh * cmesh = CreateCompMesh(data, PropertyData,porder);

  TPZGeoMesh * gmesh = cmesh->Reference();
  std::ofstream myfile("TestAnalysisTorre.txt");

  TPZNewmarkAnalysis an(cmesh,PropertyData);

  TPZSkylineStructMatrix StrMatrix(an.Mesh());
  an.SetStructuralMatrix(StrMatrix);
  TPZStepSolver<STATE> stepP;
  stepP.SetDirect( ECholesky  );
  an.SetSolver(stepP);

  an.Solution().Zero();
  an.LoadSolution();

  bool modal = true;
  REAL PeriodoVibracao = 1.;
  TPZFMatrix<STATE> AutoVetor;
  if(modal){//modal analysis

    PropertyData->SetComputeStiffness();
    an.Assemble();
    TPZAutoPointer< TPZMatrix<STATE> > K = an.Solver().Matrix()->Clone();

    PropertyData->SetComputeMass();
    an.Assemble();
    TPZAutoPointer< TPZMatrix<STATE> > M = an.Solver().Matrix()->Clone();

    TSWXPowerEigenvalueMethod powerIteration;
    TPZVec<int> Restraints(0);
    GetEquationsRestrained(*cmesh,Restraints);
    REAL EigenValue;
    TPZFMatrix<STATE> EigenVector;
    powerIteration.EigenModes(K, M, Restraints, 1e-6, 100,
                              EigenValue, EigenVector);

    an.Solution() = EigenVector;
    an.LoadSolution();
    {
      TSWXExportFrameMesh graph;
      REAL time = 0.;
      graph.GenerateGraphMesh(*cmesh,time);
      std::ofstream myfileVTK("EigenVectorTorre.vtk");
      graph.Mesh().ToParaview(myfileVTK);
      an.Solution().Zero();
      an.LoadSolution();
    }

    REAL omega = sqrt(1./EigenValue);
    REAL f = omega/(2.*M_PI);
    PeriodoVibracao = 1./f;
    AutoVetor = EigenVector;
    std::string mess;
    mess = "T = "; mess += PeriodoVibracao; mess += " s";
	std::cout << mess;
  }

  //Newmark:
  int neq = an.Mesh()->NEquations();
  TPZFMatrix<STATE> Displacement(neq,1,0.), Velocity(neq,1,0.), Acceleration(neq,1,0.);
  if(AutoVetor.Rows() == neq){
    Displacement = AutoVetor;
    Displacement *= 1e0;
  }
  an.SetInitialSolution(Displacement, Velocity, Acceleration);

  REAL alpha = 0.05;//alpha = 0, conserva energia; alpha = 0.05 é dissipaçao numerica otima
  REAL NewmarkBeta = (1.+alpha)*(1.+alpha)/4.;
  REAL NewmarkGamma = 1./2.+alpha;
  an.SetNewmarkParameters(NewmarkBeta, NewmarkGamma );

  an.TimeSteps(myfile, PeriodoVibracao/10., 0., 50, 2e-3, 60, "NewmarkTorre");

  delete cmesh;
  delete gmesh;

  return NULL;
}//void

TPZCompMesh * CreateCompMesh(TTowerData &data, TPZAutoPointer< TPZEulerBernoulliBeamData > & PropertyData,int porder){
  //problem data
  PropertyData = new TPZEulerBernoulliBeamData();

  TPZFMatrix<STATE> g(3,1,0.);
  g(2,0) = -9.8066502;
  PropertyData->SetGravity(g);

  //Material properties: MaterialId -> data
  std::map< std::string, int > MaterialName2Id;
  {//scope only
    std::map<std::string, TTowerData::TMatProp>::const_iterator w;
    int id;
    for(id = 0, w = data.fMaterialProp.begin(); w != data.fMaterialProp.end(); w++, id++){
      std::string name = w->first;
      TTowerData::TMatProp prop = w->second;
      TPZEulerBernoulliBeamData::MaterialProperties matProp;
      matProp.fE = prop.fE;
      matProp.fG = prop.fG;
      matProp.fRho = prop.fRho;
      PropertyData->fMaterialProp[id] = matProp;
      MaterialName2Id[ name ] = id;
    }//for
  }

  ///Section properties: SectionId -> data
  std::map< std::string, int > SectionName2Id;
  {//scope only
    std::map<std::string, TTowerData::TSecProp>::const_iterator w;
    int id;
    for(id = 0, w = data.fSectionProp.begin(); w != data.fSectionProp.end(); w++, id++){
      std::string name = w->first;
      TTowerData::TSecProp prop = w->second;
      TPZEulerBernoulliBeamData::SectionProperties secProp;
      secProp.fA = prop.fA;
      secProp.fIy = prop.fIy;
      secProp.fIz = prop.fIz;
      secProp.fJt = prop.fJt;
      secProp.fIp = prop.fJt;//! identical values fIp and fJt which is only true for circular cross section
      PropertyData->fSectionProp[ id ] = secProp;
      SectionName2Id[ name ] = id;
    }//for
  }

  //creating mesh
  TPZGeoMesh * gmesh = new TPZGeoMesh();

  //nodes
  std::map<int, int> JointId2PZIndex;
  {//scope only
    std::map<int, TTowerData::TJoint>::const_iterator w;
    for(w = data.fJoints.begin(); w != data.fJoints.end(); w++){
      const int nodindex = gmesh->NodeVec().AllocateNewElement();
      TPZManVector<REAL,3> nodeCoord(3);
      nodeCoord[0] = w->second.x;
      nodeCoord[1] = w->second.y;
      nodeCoord[2] = w->second.z;
      const int nodeID = w->first;
      gmesh->NodeVec()[nodindex].Initialize(nodeID,nodeCoord,*gmesh);
      JointId2PZIndex[ nodeID ] = nodindex;
    }//for
  }

  std::map<int,int> FrameID2PZIndex;
  int materialId = 3;//we do not use TPZMaterial in TPZEulerBernoulliBeam
  //elements
  {//scope only
    std::map<int, TTowerData::TIncid>::const_iterator w;
    for(w = data.fFrameIncid.begin(); w != data.fFrameIncid.end(); w++){
      const int FrameID = w->first;
      const int JointI = w->second.fJointI;
      const int JointJ = w->second.fJointJ;

    	TPZManVector<int64_t,8> incid(2);
      incid[0] = JointId2PZIndex[ JointI ];
      incid[1] = JointId2PZIndex[ JointJ ];
      int64_t index;
      gmesh->CreateGeoElement(EOned,incid,materialId,index);
      FrameID2PZIndex[ FrameID ] = index;
    }//for
  }

  //Joint restraints - supports
  int materialIdBC = -3;//we do not use TPZMaterial in TPZEulerBernoulliBeam
  std::map<int, int> JointRestrID2PZIndex;
  {//scope only
    std::map<int, TTowerData::TJointRestaint>::const_iterator w;
    for(w = data.fJointRestr.begin(); w != data.fJointRestr.end(); w++){
      int JointID = w->first;
      TPZManVector<int64_t,8> incid(1);
      incid[0] = JointId2PZIndex[ JointID ];
      int64_t index;
      gmesh->CreateGeoElement(EPoint,incid,materialIdBC,index);
      JointRestrID2PZIndex[ JointID ] = index;
    }
  }

	gmesh->BuildConnectivity();

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

  //elements
  {//scope only
    std::map<int, TTowerData::TIncid>::const_iterator w;
    for(w = data.fFrameIncid.begin(); w != data.fFrameIncid.end(); w++){
      const int FrameID = w->first;
      const int gelIndex = FrameID2PZIndex[ FrameID ];
      TPZGeoEl * gel = gmesh->ElementVec()[ gelIndex ];
      if(!gel) DebugStop();
      int64_t index;
      TPZEulerBernoulliBeam * cel = new TPZEulerBernoulliBeam(*cmesh, gel, index);

      //section
      std::string sectionName = data.fFrameSection[ FrameID ];
      const int SectionId = SectionName2Id[ sectionName ];

      //material
      TTowerData::TSecProp sectionProp =  data.fSectionProp[sectionName];
      std::string materialName = sectionProp.fMaterial;
      const int MaterialId = MaterialName2Id[ materialName ];

      //section rotation - if it exists
      REAL alfa;
      std::map<int, double>::const_iterator it = data.fFrameSectionRotation.find( FrameID );
      if(it == data.fFrameSectionRotation.end()){
        alfa = 0;
      }
      else{
        alfa = it->second * M_PI/180.;//converting from º do rad
      }

      REAL fabStrainError = 0.;
      cel->SetPropertyData( PropertyData, MaterialId, SectionId, alfa, fabStrainError );
    }
  }

  //Joint restraints - supports
  {//scope only
    std::map<int, TTowerData::TJointRestaint>::const_iterator w;
    for(w = data.fJointRestr.begin(); w != data.fJointRestr.end(); w++){
      int JointID = w->first;
      int gelindex = JointRestrID2PZIndex[ JointID ];
      TPZGeoEl * gel = gmesh->ElementVec()[gelindex];
      if(!gel) DebugStop();
      int64_t index;
      TPZEulerBernoulliBC * cel = new TPZEulerBernoulliBC(*cmesh, gel, index);
      cel->SetData( PropertyData );
      TTowerData::TJointRestaint jointRestData = w->second;
      TPZVec< TPZEulerBernoulliBC::BCVal > bcdata(6);
      if(  jointRestData.fU1 ) bcdata[0] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      if(  jointRestData.fU2 ) bcdata[1] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      if(  jointRestData.fU3 ) bcdata[2] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      if(  jointRestData.fR1 ) bcdata[3] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      if(  jointRestData.fR2 ) bcdata[4] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      if(  jointRestData.fR3 ) bcdata[5] = TPZEulerBernoulliBC::BCVal( TPZEulerBernoulliBC::ESupport, 0. );
      cel->SetBCVal(bcdata);
    }
  }

  // Introduzindo um material apenas para o Calculo do Residuo funcionar
  TPZVec<REAL> nullforce(3, 0.);
  TPZElasticity3DGD * elastMat = new TPZElasticity3DGD(materialId, 33488., 0.2, nullforce, 0., 0., 0.);
  cmesh->InsertMaterialObject(elastMat);
  // Tambem material BC
  TPZFMatrix<STATE> val1(3, 3, 0.), val2(3, 1, 0.);
  TPZMaterial* bcmat = elastMat->CreateBC(elastMat, materialIdBC, 0, val1, val2);
  cmesh->InsertMaterialObject(bcmat);

	cmesh->SetDefaultOrder( porder );
  cmesh->SetDimModel(3);

	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();

	return cmesh;
}//CreateCompMesh

void GetEquationsRestrained(TPZCompMesh & cmesh,TPZVec<int> &Restraints){
  std::set<int> eqs;
  TPZManVector<int> indices;
  for(int iel = 0; iel < cmesh.NElements(); iel++){
    TPZEulerBernoulliBC * bcel = dynamic_cast< TPZEulerBernoulliBC* >(cmesh.ElementVec()[iel]);
    if(!bcel) continue;
    bcel->GetEquationIndices(indices);
    for(int i = 0; i < 6; i++){
      if( bcel->GetBCVal(i).fType == TPZEulerBernoulliBC::ESupport){
        eqs.insert(indices[i]);
      }
    }
  }
  Restraints.Resize(eqs.size());
  std::set<int>::iterator w;
  int i;
  for(i = 0, w = eqs.begin(); w != eqs.end(); w++, i++){
    Restraints[i] = *w;
  }
}//void


bool TSWXPowerEigenvalueMethod::EigenModes(TPZAutoPointer<TPZMatrix<STATE> > Stiffness,
	TPZAutoPointer<TPZMatrix<STATE> > &Mass,
                  const TPZVec<int> & Restraints, REAL tol, int niter,
				  REAL & EigenValue, TPZFMatrix<STATE> & EigenVector){

  const int neq = Mass->Rows();
  //criando vetores auxiliares
  TPZFMatrix<STATE> u(neq, 1);
  EigenVector.Redim(neq,1);

  std::list<int64_t> singular;
  Stiffness->Decompose_Cholesky(singular);
  if(singular.size()){
    std::cout << ( "Stiffness->Decompose_Cholesky has found singular modes" );
  }

  //chute inicial - valores aleatórios
  for(int i = 0; i < neq; i++){
    EigenVector(i,0) = static_cast<double>( rand() % 100 )/100.;
  }

  int iter = 0;
  REAL lambda = 0.0, diff = 2.0 * tol + 1.;

  while (diff > tol && iter < niter) {
    iter++;
    Mass->Multiply(EigenVector, u);
    //Boundary restriction
    for(int i = 0; i < Restraints.NElements(); i++){
      int eq = Restraints[i];
      u( eq, 0 ) = 0.;           //talvez nao precise disto aqui por causa do bignumber nas condicoes dirichlet
    }
    Stiffness->Solve_Cholesky(&u);
    lambda = Norm(u);
    u *= 1.0/lambda;
    EigenVector -= u;
    diff = Norm(EigenVector);
    EigenVector = u;
  }

  EigenValue = lambda;

#ifdef DEBUG
  std::string mess;
  mess = "EigenValue = ";
  mess += EigenValue;
  mess += " , diff = ";
  mess += diff;
#endif

  if(tol > diff) return false;
  else return true;

}
