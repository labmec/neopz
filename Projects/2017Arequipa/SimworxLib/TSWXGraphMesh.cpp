
#include "TSWXGraphMesh.h"
#include "pzerror.h"
#include "fstream"

/////////////////// TSWXGraphSol //////////////////////

int TSWXGraphSol::Dimension() const{
	if(fData.size()) return fData[0].size();
	else return 0;
}

unsigned int TSWXGraphSol::Size() const
{
	return fData.size();
}


void TSWXGraphSol::Read(std::istream & file){
  int soltype;
  file >> soltype;
  fSolType = (ESWXVtkSolType)(soltype);
  file >> fTitle;
  int n;
  file >> n;
  fData.resize(n);
  for(unsigned int i = 0; i < fData.size(); i++){
//    fData[i].Read(file);
  }
}

void TSWXGraphSol::Write(std::ostream &file) const{
  file << (int)(fSolType) << "\n";
  file << fTitle << "\n";
  file << fData.size() << "\n";
  for(unsigned int i = 0; i < fData.size(); i++){
//    fData[i].Write(file);
  }
  file << "\n";
}

/////////////// TSWXGraphMeshSol //////////////////////

TSWXGraphMeshSol::TSWXGraphMeshSol() : fAllSol() {

}

TSWXGraphMeshSol::~TSWXGraphMeshSol(){

}

void TSWXGraphMeshSol::AddSol(double time, const TPZVec< TSWXGraphSol >  &sol){
  this->fAllSol[ time ] = sol;
}

unsigned int TSWXGraphMeshSol::NTimeSteps() const{
  return fAllSol.size();
}

void TSWXGraphMeshSol::GetSolData(int iStep, double &time, TPZVec< TSWXGraphSol >  &sol) const{
  sol = GetSolData(iStep,time);
}

const TPZVec< TSWXGraphSol >  & TSWXGraphMeshSol::GetSolData(int iStep, double &time) const{
  if(iStep < 0 || iStep >= this->NTimeSteps()){
    DebugStop();
  }

  std::map< double, TPZVec< TSWXGraphSol > >::const_iterator w = this->fAllSol.begin();
  for(int i = 0; i < iStep; i++) w++; ///Tiago: acho que dá w += iStep, mas nao tenho certeza
  time = w->first;
  return w->second;
}

void TSWXGraphMeshSol::Clear(){
  this->fAllSol.clear();
}

void TSWXGraphMeshSol::Read(std::istream & file){
  std::string aux;
  file >> aux;

  this->Clear();
  int n;
  file >> n;
  for(int i = 0; i < n; i++){
    double time;
    file >> time;
    TPZVec< TSWXGraphSol > sol;
//    ReadVec(file, sol);
    fAllSol[time] = sol;
  }
}

void TSWXGraphMeshSol::Write(std::ostream &file) const{
	file << "TSWXGraphMeshSol::Write\n";
	file << fAllSol.size() << "\n";
	std::map< double, TPZVec< TSWXGraphSol > >::const_iterator w;
	for(w = fAllSol.begin(); w != fAllSol.end(); w++){
		file << w->first << "\n";
//		WriteVec( file, w->second );
	}
}
TPZVec<std::string> TSWXGraphMeshSol::GetSolutionsNames() const
{
	TPZVec<std::string> names(fAllSol.begin()->second.size());
	for(int i = 0; i < fAllSol.begin()->second.size(); i++)
	{
		names[i] = fAllSol.begin()->second[i].fTitle;
  }
	return names;
}

void TSWXGraphMeshSol::
				GetSolution(std::string SolName, std::map<double , TSWXGraphSol> & solMap) const
{
	std::map< double, TPZVec< TSWXGraphSol > >::const_iterator it;
	int solId = -1;
	for(int i = 0; i < fAllSol.begin()->second.size(); i++)
	{
		solId = i;
		if(SolName == fAllSol.begin()->second[i].fTitle) break;
	}
	for(it = fAllSol.begin(); it != fAllSol.end(); it++)
	{
		double instant = it->first;
		TSWXGraphSol sol = it->second[solId];
		//std::pair<double,  TSWXGraphSol > item(instant, sol);
		solMap[instant] = sol;
	}
}


/////////////////  TSWXGraphEl  //////////////////

void TSWXGraphEl::Read(std::istream & file){
//  fIncid.Read(file);
	int eltype;
  file >> eltype;
  fElType = (EGraphElType)(eltype);
  file >> fMatId;
}

void TSWXGraphEl::Write(std::ostream &file) const{
//  fIncid.Write(file);
  file << (int)(fElType) << "\n";
  file << fMatId << "\n";
}

//////////////   TSWXGraphMesh  /////////////////

void TSWXGraphMesh::Read(std::istream & file){

  unsigned int n;

  file >> n;
  fNodes.resize(n);
  for(unsigned int i = 0; i < n; i++){
//    fNodes[i].Read(file);
  }

  file >> n;
  fElem.resize(n);
  for(unsigned int i = 0; i < n; i++){
    fElem[i].Read(file);
  }

  fSol.Read(file);

  file >> n;
  fMaterialLabels.clear();
  for(unsigned int i = 0; i < n; i++){
    int id;
    std::string label;
    file >> id >> label;
    fMaterialLabels[id] = label;
  }

}

void TSWXGraphMesh::Write(std::ostream &file) const{

  file << fNodes.size() << "\n";
  for(unsigned int i = 0; i < fNodes.size(); i++){
 //   fNodes[i].Write(file);
  }

  file << fElem.size() << "\n";
  for(unsigned int i = 0; i < fElem.size(); i++){
    fElem[i].Write(file);
  }

  fSol.Write(file);

  file << fMaterialLabels.size() << "\n";
  std::map< int, std::string >::const_iterator w;
  for( w = fMaterialLabels.begin(); w != fMaterialLabels.end(); w++){
    file << w->first << "\t" << w->second << "\n";
  }

}///void

bool TSWXGraphMesh::HasNodalSolution() const{
  if(fSol.NTimeSteps()){
    double time;
    const TPZVec< TSWXGraphSol > & sol = fSol.GetSolData(0,time);
    const int n = sol.size();
    for(int i = 0; i < n; i++){
      if(sol[i].fSolType == ENodeSolution) return true;
    }
    return false;
  }
  else return false;
}

bool TSWXGraphMesh::HasCellSolution() const{
  if(fSol.NTimeSteps()){
    double time;
    const TPZVec< TSWXGraphSol > & sol = fSol.GetSolData(0,time);
    const int n = sol.size();
    for(int i = 0; i < n; i++){
      if(sol[i].fSolType == ECellSolution) return true;
    }
    return false;
  }
  else return false;
}

int TSWXGraphMesh::VTKCellType( EGraphElType enumType ) const{
  switch(enumType) {
    case EvtkQuadraticEdge:
      return 21;
      break;
    case EvtkBiQuadraticTriangle:
      return 34;
      break;
    case EvtkBiQuadraticQuad:
      return 28;
      break;
    case EvtkQuadraticTetra:
      return 24;
      break;
    case EvtkQuadraticPyramid:
      return 27;
      break;
    case EvtkQuadraticWedge:
      return 26;
      break;
    case EvtkQuadraticHexahedron:
      return 25;
      break;
    case EvtkLinearHexahedron:
      return 12;
      break;
    default:
      DebugStop();
      return -1;
  }
  DebugStop();
  return -1;
}

void TSWXGraphMesh::ToParaview(std::ostream &file, int mathStyleIStep) const{

  int iStep;
  if(mathStyleIStep < 0){
    iStep = fSol.NTimeSteps() + mathStyleIStep;
  }
  else{
    iStep = mathStyleIStep;
  }

	//Header
	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "Simworx Eng P&D" << std::endl;
	file << "ASCII" << std::endl << std::endl;

	file << "DATASET UNSTRUCTURED_GRID" << std::endl;
	file << "POINTS ";

  const int nnodes = this->fNodes.size();
	file << nnodes << " float" << std::endl;
  for(int i = 0; i < nnodes; i++){
    int n = this->fNodes[i].size();
    for(int j = 0; j < n; j++){
      file << this->fNodes[i][j] << " ";
    }
    file << "\n";
  }///nodes

  const int nelem = this->fElem.size();
  int nvalues = 0;
  for(int i = 0; i < nelem; i++){
    nvalues += this->fElem[i].fIncid.size() + 1;
  }

	file << "CELLS " << nelem << " ";
	file << nvalues << std::endl;
  for(int i = 0; i < nelem; i++){
    int n = this->fElem[i].fIncid.size();
    file << n << " ";
    for(int j = 0; j < n; j++){
      file << this->fElem[i].fIncid[j] << " ";
    }
    file << "\n";
  }///elements

  file << "CELL_TYPES " << nelem << "\n";
  for(int i = 0; i < nelem; i++){
    file << this->VTKCellType( this->fElem[i].fElType ) << "\n";
  }

  ///malha vazia
  if(nelem == 0) return;

  if(this->HasNodalSolution()){
    file << "POINT_DATA " << nnodes << "\n";
    double time;
    const TPZVec< TSWXGraphSol > & sol = fSol.GetSolData(iStep,time);
    const int nSol = sol.size();
    for(int i = 0; i < nSol; i++){
      if(sol[i].fSolType != ENodeSolution) continue;
      if(sol[i].Dimension() == 1){///escalar
        file << "SCALARS " << sol[i].fTitle << " float\n";
        file << "LOOKUP_TABLE default\n";
        for(unsigned int j = 0; j < sol[i].fData.size(); j++){
          file << sol[i].fData[j][0] << "\n";
        }///j
      }///fim escalar
      else{///vetorial
        if(sol[i].Dimension() != 3){
          {
            std::ofstream myfile("c:\\Temp\\TSWXGraphMeshErro.txt");
            sol[i].Write(myfile);
          }
          DebugStop();
        }
        if(sol[i].fSolType != ENodeSolution) continue;
        file << "VECTORS " << sol[i].fTitle << " float\n";
        for(unsigned int j = 0; j < sol[i].fData.size(); j++){
          file << sol[i].fData[j][0] << "\t"
               << sol[i].fData[j][1] << "\t"
               << sol[i].fData[j][2] << "\n";
        }///j
      }///fim vetorial

    }///i solutions
  }///has nodal solution

  if(this->HasCellSolution()){
    file << "CELL_DATA " << nelem << "\n";
    double time;
    const TPZVec< TSWXGraphSol > & sol = fSol.GetSolData(iStep,time);
    const int nSol = sol.size();
    for(int i = 0; i < nSol; i++){
      if(sol[i].fSolType != ECellSolution) continue;
      if(sol[i].Dimension() == 1){///escalar
        file << "SCALARS " << sol[i].fTitle << " float\n";
        file << "LOOKUP_TABLE default\n";
        for(unsigned int j = 0; j < sol[i].fData.size(); j++){
          file << sol[i].fData[j][0] << "\n";
        }///j
      }///fim escalar
      else{///vetorial
        if(sol[i].Dimension() != 3) DebugStop();
        if(sol[i].fSolType != ECellSolution) continue;
        file << "VECTORS " << sol[i].fTitle << " float\n";
        for(unsigned int j = 0; j < sol[i].fData.size(); j++){
          file << sol[i].fData[j][0] << "\t"
               << sol[i].fData[j][1] << "\t"
               << sol[i].fData[j][2] << "\n";
        }///j
      }///fim vetorial

    }///i solutions
  }

}///void


/*_di_IXMLNode TSWXGraphMesh::ReadMe(_di_IXMLNode & fatherNode, const UnicodeString & nodeName)
{
  _di_IXMLNode myNode = fatherNode->ChildNodes->FindNode(nodeName);
  if(myNode)
  {
    std::istringstream vtkStr;
    std::string conv(swx::wstring2string(myNode->Text.c_str()));
    vtkStr.str(conv);

    this->Read(vtkStr);
  }
  else
  {
    this->Reset();
  }

  return myNode;
}

_di_IXMLNode TSWXGraphMesh::WriteMe(_di_IXMLNode & fatherNode, const UnicodeString & nodeName) const
{
  std::stringstream vtkStr;
  this->Write(vtkStr);

  _di_IXMLNode vtkNode = fatherNode->AddChild(nodeName);
  vtkNode->Text = vtkStr.str().c_str();

  return vtkNode;
}

*/

