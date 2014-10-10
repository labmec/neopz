/**
 * @file
 * @brief Contains the implementation of the TPZGeoElRefLess methods.
 */

#ifndef PZGEOELREFLESS_H_H
#define PZGEOELREFLESS_H_H

#include "pzgeoelrefless.h"

#include <sstream>

#include "pzlog.h" // test 
#ifdef LOG4CXX
static LoggerPtr loggerrefless(Logger::getLogger("pz.mesh.tpzgeoelrefless"));
#endif

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess():TPZGeoEl(){
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSideIndex();
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(const TPZGeoElRefLess<TGeo>  &gel):TPZGeoEl(gel), fGeo(gel.fGeo){
	int i;
	for(i=0;i<TGeo::NSides;i++){
		TPZGeoElSide thisside(this->fNeighbours[i], this->Mesh());
		TPZGeoElSide gelside(gel.fNeighbours[i], this->Mesh());
		thisside.SetConnectivity(gelside);
		//     fNeighbours[i].SetConnectivity(gel.fNeighbours[i]);
	}
	//  fSubElement = -1;
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::~TPZGeoElRefLess(){
	//RemoveConnectivities();
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZVec<long> &nodeindices,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(matind,mesh), fGeo(nodeindices) {
	
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSideIndex();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TGeo &geo,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(matind,mesh), fGeo(geo) {
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSideIndex();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZVec<long> &nodeindices,int matind,TPZGeoMesh &mesh, long &index) :
TPZGeoEl(matind,mesh,index) , fGeo(nodeindices) 
{
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSideIndex();
    fGeo.Initialize(this);
}

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(long id,TPZVec<long> &nodeindexes,int matind,TPZGeoMesh &mesh) :
TPZGeoEl(id,matind,mesh) , fGeo(nodeindexes) {
	int i;
	for(i=0;i<TGeo::NSides;i++)fNeighbours[i] = TPZGeoElSideIndex();
    fGeo.Initialize(this);
}

template<class TGeo>
long
TPZGeoElRefLess<TGeo>::NodeIndex(int node) const {
	if(node<0 || node>=fGeo.NNodes) return -1;
	return fGeo.fNodeIndexes[node];
}

template<class TGeo>
long
TPZGeoElRefLess<TGeo>::SideNodeIndex(int side,int node) const {
	if(side<0 || side>(TGeo::NSides - 1) || node<0) {
		PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
		return -1;
	}
	return fGeo.fNodeIndexes[TGeo::SideNodeLocId(side,node)];
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideNodeLocIndex(int side,int node) const{
	
	if(side<0 || side>(TGeo::NSides - 1) || node<0) {
		PZError << "TPZGeoElRefLess::SideNodeIndex. Bad parameter side.\n";
		return -1;
	}
	return TGeo::SideNodeLocId(side,node);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::SetSubElement(int id, TPZGeoEl *el){
	if (id!=0)
		PZError << "TPZGeoElRefLess<TGeo>::SetSubElement - Fodeu!\n";
	else PZError << "TPZGeoElRefLess<TGeo>::SetSubElement - Por enquanto eu no faco nada!\n";
	return;
}

template<class TGeo>
TPZIntPoints *TPZGeoElRefLess<TGeo>::CreateSideIntegrationRule(int side, int order){
	return TGeo::CreateSideIntegrationRule(side,order);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NNodes() const {
	return TGeo::NNodes;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NCornerNodes() const{
	return TGeo::NCornerNodes;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSides() const{
	return TGeo::NSides;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideNodeLocId(int side, int node) const {
	return TGeo::SideNodeLocId(side,node);
}

template<class TGeo>
REAL
TPZGeoElRefLess<TGeo>::RefElVolume(){
	return TGeo::RefElVolume();
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSideNodes(int side) const{
	return TGeo::NSideNodes(side);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::MidSideNodeIndex(int side,long &index) const{
	//TRef::MidSideNodeIndex(this,side,index);
	index = -1;
	if(side<0 || side>NSides()-1) {
		PZError << "TPZGeoElRefLess<TGeo>::MidSideNodeIndex. Bad parameter side = " << side << std::endl;
		return;
	}
	if(side<NNodes()) {//o n�medio do lado 0 �o 0 etc.
		index = NodeIndex(side);
		return;
	}
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideIsUndefined(int side){
	if (side < 0 || side > NSides()){
		PZError << "TPZGeoElRefLess<TGeo>::SideIsUndefined - bad side: " << side << std::endl;
	}
	return (fNeighbours[side].Side() == -1);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSubElements() const {
	//return TRef::NSubEl;
	return 0;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::NSideSubElements(int side) const {
	return 0;
}


template<class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TGeo>::CreateBCGeoEl(int side, int bc){
	TPZGeoEl * result = fGeo.CreateBCGeoEl(this,side,bc);
	//result->BuildBlendConnectivity();
	return result;  
}

template<class TGeo>
TPZGeoEl * TPZGeoElRefLess<TGeo>::CreateGeoElement(MElementType type,
												   TPZVec<long>& nodeindexes,
												   int matid,
												   long& index)
{
	return fGeo.CreateGeoElement(*Mesh(),type,nodeindexes,matid,index);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::SetNodeIndex(int i,long nodeindex){
	if(i<0 || i>(TGeo::NNodes - 1)){
		std::cout << "TPZGeoElRefLess::SetNodeIndex index error i = " << i << std::endl;
		return;
	}
	fGeo.fNodeIndexes[i] = nodeindex;
}

template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::SideToSideTransform(int sidefrom,int sideto){
	return TGeo::SideToSideTransform(sidefrom,sideto);
}

template<class TGeo>
TPZGeoEl *
TPZGeoElRefLess<TGeo>::SubElement(int is) const {
	if(is<0 || is>1){//(TRef::NSubEl - 1)){
		std::cout << "TPZGeoElRefLess::SubElement index error is= " << is << std::endl;
	}
	//  return fSubEl[is];
	return 0;
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::SideDimension(int side) const {
	return TGeo::SideDimension(side);
}

template<class TGeo>
int
TPZGeoElRefLess<TGeo>::Dimension() const {
	return TGeo::Dimension;
}

template<class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TGeo>::HigherDimensionSides(int side,int targetdimension){
	std::cout << "TPZGeoElRefLess::HigherDimensionSides nao deve ser usado\n";
	return TPZGeoElSide();
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::AllHigherDimensionSides(int side,int targetdimension,TPZStack<TPZGeoElSide> &elsides){
	TPZStack<int> highsides;
	TGeo::HigherDimensionSides(side,highsides);
	int i,size = highsides.NElements();
	for(i=0;i<size;i++) {
		if(SideDimension(highsides[i]) == targetdimension) {
			elsides.Push(TPZGeoElSide(this,highsides[i]));
		}
	}
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::LowerDimensionSides(int side,TPZStack<int> &smallsides){
	TGeo::LowerDimensionSides(side,smallsides);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::BuildTransform(int side, TPZGeoEl *father,TPZTransform &t){
	BuildTransform2(side,father,t);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::Jacobian(TPZVec<REAL> &coordinate,TPZFMatrix<REAL> &jac,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const {
	fGeo.Jacobian(*this,coordinate,jac,axes,detjac,jacinv);
	
    
#ifdef DEBUG
    TPZManVector<REAL,3> minx(3,0.),maxx(3,0.);
    int ncorners = NNodes();
    TPZFNMatrix<54,REAL> cornerco(3,ncorners);
    fGeo.CornerCoordinates(*this,cornerco);
    long spacedim = cornerco.Rows();
    
    for (long j=0; j<spacedim; j++) {
        minx[j] = cornerco(j,0);
        maxx[j] = cornerco(j,0);
    }
    for (int i=0; i<ncorners; i++) {
        for (int d=0; d<spacedim; d++) {
            minx[d] = minx[d] < cornerco(d,i) ? minx[d] : cornerco(d,i);
            maxx[d] = maxx[d] > cornerco(d,i) ? maxx[d] : cornerco(d,i);
        }
    }
    REAL delx=0.;
    for (int i=0; i<ncorners; i++) {
        for (long d=0; d<spacedim; d++) {
            delx = delx < maxx[d]-minx[d] ? maxx[d]-minx[d] : delx;
        }
    }

	if(TGeo::Dimension != 0 && (IsZero(delx) || IsZero(detjac/(delx*delx)))){
		std::stringstream sout;
		sout << "Jacobiano nulo\n";
		LOGPZ_ERROR(loggerrefless,sout.str())
		detjac = ZeroTolerance();
	}
#endif
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::X(TPZVec<REAL> &coordinate,TPZVec<REAL> &result) const {
	fGeo.X(*this,coordinate,result);
}

template<class TGeo>
bool TPZGeoElRefLess<TGeo>::IsLinearMapping() const 
{ 
	return fGeo.IsLinearMapping();
}

template<class TGeo>
bool TPZGeoElRefLess<TGeo>::IsGeoBlendEl() const 
{ 
	return fGeo.IsGeoBlendEl();
}

template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::BuildTransform2(int side, TPZGeoEl * father, TPZTransform &t)
{
	if(this == father) return t;
	TPZGeoEl *myfather = Father();
	if(side<0 || side>(TGeo::NSides-1) || !myfather){
		PZError << "TPZGeoElRefLess::BuildTransform2 side out of range or father null\n";
		return TPZTransform(0,0);
	}
	TPZGeoElSide fathloc = Father2(side);
	int son = WhichSubel();
	TPZTransform trans=myfather->GetTransform(side,son);
	trans = trans.Multiply(t);
	if(fathloc.Element() == father) return trans;
	trans = myfather->BuildTransform2(fathloc.Side(),father,trans);
	return trans;
}

template<class TGeo>
TPZTransform
TPZGeoElRefLess<TGeo>::GetTransform(int /*side*/,int /*son*/){
    PZError << "TPZGeoElRefLess<TGeo>::GetTransform::Never should be called\n";
    return TPZTransform(0,0);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::CenterPoint(int side, TPZVec<REAL> &cent) const {
	TGeo::CenterPoint(side,cent);
}

template<class TGeo>
TPZGeoElSide
TPZGeoElRefLess<TGeo>::Father2(int side)
{
	//std::cout << " Father2 teste Cedric: 08/05/2003\n";
	TPZGeoEl *father = Father();
	if(!father) return TPZGeoElSide();
	int son = WhichSubel();
	if(son<0) return TPZGeoElSide();
	int fathsid = father->FatherSide(side,son);
	return TPZGeoElSide(father,fathsid);
}

template<class TGeo>
void
TPZGeoElRefLess<TGeo>::GetSubElements2(int side, TPZStack<TPZGeoElSide> &subel) const{
	return;
}

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Read(TPZStream &buf, void *context){
	TPZGeoEl::Read(buf,context);
    fGeo.Read(buf,context);
#ifdef DEBUG
    long NNodes = Mesh()->NodeVec().NElements();
    for (int i=0; i< TGeo::NNodes; i++) {
        if (fGeo.fNodeIndexes[i]<0 || fGeo.fNodeIndexes[i] > NNodes) {
            DebugStop();
        }
    }
#endif
    
	int i, n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i].Read(buf);
#ifdef DEBUG
        int nel = Mesh()->NElements();
        if (this->fNeighbours[i].ElementIndex() < 0 || this->fNeighbours[i].Side() < 0 || this->fNeighbours[i].Side() >=27
            || this->fNeighbours[i].ElementIndex() >= nel) {
            DebugStop();
        }
#endif
	}
}//Read

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Write(TPZStream &buf, int withclassid){
	TPZGeoEl::Write(buf,withclassid);
    fGeo.Write(buf);
	int i, n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i].Write(buf);
	}
}//Write

template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess(TPZGeoMesh &DestMesh, const TPZGeoElRefLess &cp):TPZGeoEl(DestMesh, cp), fGeo(cp.fGeo) {
	int i;
	const int n = TGeo::NSides;
	for(i = 0; i < n; i++){
		this->fNeighbours[i] = cp.fNeighbours[i];
	}
}


template<class TGeo>
TPZGeoElRefLess<TGeo>::TPZGeoElRefLess( TPZGeoMesh &DestMesh,
									   const TPZGeoElRefLess &cp,
									   std::map<long,long> & gl2lcNdMap,
									   std::map<long,long> & gl2lcElMap ) :
TPZGeoEl(DestMesh, cp, gl2lcElMap), fGeo(cp.fGeo, gl2lcNdMap)
{
	int i;
	const int n = TGeo::NSides;
	
	for(i = 0; i < n; i++)
	{
		TPZGeoElSide neigh (cp.fNeighbours[i],cp.Mesh());
		long neighIdx = neigh.Element()->Index();
		int side = neigh.Side();
		
		while (gl2lcElMap.find(neighIdx)==gl2lcElMap.end())
		{
			neigh = neigh.Neighbour();
			neighIdx = neigh.Element()->Index();
			side = neigh.Side();
		}
		this->fNeighbours[i] = TPZGeoElSideIndex ( gl2lcElMap [ neighIdx ] , side );
	}
}

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Directions(int side, TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
{
    TPZFNMatrix<9,REAL> jac(TGeo::Dimension,TGeo::Dimension), jacinv(TGeo::Dimension,TGeo::Dimension), axes(TGeo::Dimension,3), gradx(3,TGeo::Dimension,0.);
    REAL detjac;

    this->Jacobian(pt,jac,axes,detjac,jacinv);
    
    // ou eh isso?   grad =  (jac  * axes)ˆT 
    TPZFNMatrix<9> gradxt(TGeo::Dimension,3,0.);
    for (int il=0; il<TGeo::Dimension; il++)
    {
        for (int jc=0; jc<3; jc++)
        {
            for (int i = 0 ; i<TGeo::Dimension; i++)
            {
                gradx(jc,il) += jac(i,il) * axes(i,jc);
            }
        }
    }
//    gradxt.Transpose(&gradx);
    TGeo::ComputeDirections(side, gradx, directions, sidevectors);
    
//    TPZStack<int> lowdim;
//	LowerDimensionSides(side,lowdim);
//	lowdim.Push(side);
//
//    TGeo::GetSideHDivPermutation(side);
    
}

template<class TGeo>
void TPZGeoElRefLess<TGeo>::Directions(TPZVec<REAL> &pt, TPZFMatrix<REAL> &directions)
{
    TPZFNMatrix<9,REAL> jac(TGeo::Dimension,TGeo::Dimension), jacinv(TGeo::Dimension,TGeo::Dimension), axes(TGeo::Dimension,3), gradx(3,TGeo::Dimension,0.);
    REAL detjac;
    
    this->Jacobian(pt,jac,axes,detjac,jacinv);
    
    // ou eh isso?   grad =  (jac  * axes)ˆT
    TPZFNMatrix<9> gradxt(TGeo::Dimension,3,0.);
    for (int il=0; il<TGeo::Dimension; il++)
    {
        for (int jc=0; jc<3; jc++)
        {
            for (int i = 0 ; i<TGeo::Dimension; i++)
            {
                gradx(jc,il) += jac(i,il) * axes(i,jc);
            }
        }
    }
    //    gradxt.Transpose(&gradx);
    TGeo::ComputeDirections(gradx, detjac, directions);
    
}


#include "pzgeoquad.h"

/** Compute the permutation for an HDiv side */
/*
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoQuad>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	if(side < 4 || side > 7)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
        std::cout << sout.str() << std::endl;
	}
	permutegather.Resize(3);
	long id1 = NodePtr(SideNodeLocIndex(side,0))->Id();
	long id2 = NodePtr(SideNodeLocIndex(side,1))->Id();
	if(id1<id2)
	{
		permutegather[0] = 0;
		permutegather[1] = 1;
		permutegather[2] = 2;
	}
	else
	{
		permutegather[0] = 1;
		permutegather[1] = 0;
		permutegather[2] = 2;
	}
}
*/
#include "pzgeotriangle.h"

/** Compute the permutation for an HDiv side */
/*
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoTriangle>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	if(side < 3 || side > 5)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
        std::cout << sout.str() << std::endl;
	}
	permutegather.Resize(3);
	long id1 = NodePtr(SideNodeLocIndex(side,0))->Id();
	long id2 = NodePtr(SideNodeLocIndex(side,1))->Id();
	if(id1<id2)
	{
		permutegather[0] = 0;
		permutegather[1] = 1;
		permutegather[2] = 2;
	}
	else
	{
		permutegather[0] = 1;
		permutegather[1] = 0;
		permutegather[2] = 2;
	}
}
*/

/** Compute the permutation for an HDiv side */
template<class TGeo>
inline void TPZGeoElRefLess<TGeo>::HDivPermutation(int side, TPZVec<int> &permutegather)
{
	int dimension = TGeo::Dimension;
	int sidedimension = TGeo::SideDimension(side);
	
	if(dimension != sidedimension+1)
	{
		std::stringstream sout;
		sout << "HDivPermutation called with wrong side parameter " << side;
#ifdef LOG4CXX
		LOGPZ_ERROR(loggerrefless,sout.str())
#endif
	}
    
    // Douglas -- teste em 2014 09 04
    // conta o numero de lados da face
    const long nsidenodes = TGeo::NSideNodes(side);
    TPZManVector<long,4> id(nsidenodes);  // 
    
	for(int inode=0; inode<nsidenodes; inode++)
    {
        // esta parte pega os indices locais dos nos apenas da face em questao
        long nodeindex = SideNodeLocId(side, inode);
        // com base nestes indices locais, pegamos os indices globais para determinar a permutacao
        id[inode] = NodePtr(nodeindex)->Id();
    }
    
    // Esse bloco parece pegar todo os vertices do cubo para fazer a permutacao, deveria ser da face
//    TPZManVector<long,TGeo::NCornerNodes> id(TGeo::NCornerNodes);
//	for(int i=0; i<TGeo::NCornerNodes; i++)
//    {
//        long nodeindex = fGeo.fNodeIndexes[i];
//        id[i] = Mesh()->NodeVec()[nodeindex].Id();
//    }
    
    MElementType sidetype = TGeo::Type(side);
    int transformid;
    switch (sidetype) {
        case EOned:
            transformid = pztopology::TPZLine::GetTransformId(id);
            pztopology::TPZLine::GetSideHDivPermutation(transformid, permutegather);
            break;
        case EQuadrilateral:
            transformid = pztopology::TPZQuadrilateral::GetTransformId(id);
            pztopology::TPZQuadrilateral::GetSideHDivPermutation(transformid, permutegather);
            break;
        case ETriangle:
            transformid = pztopology::TPZTriangle::GetTransformId(id);
            pztopology::TPZTriangle::GetSideHDivPermutation(transformid, permutegather);
            break;
        default:
            DebugStop();
            break;
    }
#ifdef LOG4CXX
    if (loggerrefless->isDebugEnabled()) {
        std::stringstream sout;
        sout << "side = " << side << " transform id " << transformid << " permutegather " << permutegather;
        LOGPZ_DEBUG(loggerrefless, sout.str())
    }
#endif
}

//HDiv
template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoQuad>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
    fGeo.VecHdiv(*this,normalvec,sidevector);
}

template<>
inline void TPZGeoElRefLess<pzgeom::TPZGeoTriangle>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
	fGeo.VecHdiv(*this,normalvec,sidevector);
}


template<class TGeo>
inline void TPZGeoElRefLess<TGeo>::VecHdiv(TPZFMatrix<REAL> &normalvec,TPZVec<int> &sidevector )
{
    PZError << __PRETTY_FUNCTION__ << " nao esta implementado\n";
}

#endif
