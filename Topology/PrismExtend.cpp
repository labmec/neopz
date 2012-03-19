/**
 * @file
 * @brief Contains the implementation of the Pr methods. 
 */
#include "PrismExtend.h"
#include "tpzpoint.h"

using namespace std;

#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr logger(Logger::getLogger("pz.topology.prismextend"));
#endif

namespace pztopology {
	
	template<class TFather>
	Pr<TFather>::Pr()
	{
	}
	
	template<class TFather>
	Pr<TFather>::~Pr()
	{
	}
	
	template<class TFather>
	int Pr<TFather>::NSideNodes(int side)
	{
		int ftns = side/TFather::NSides;
		if(ftns < 2) return TFather::NSideNodes(side % TFather::NSides);
		return TFather::NSideNodes(side % TFather::NSides) *2;
	}
	
	template<class TFather>
	int Pr<TFather>::SideNodeLocId(int side, int node)
	{
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		if(ftns == 0) return TFather::SideNodeLocId(side,node);
		if(ftns == 1) return TFather::SideNodeLocId(side-TFather::NSides,node)+TFather::NCornerNodes;
		int fatherlevel = node/TFather::NSideNodes(fatherside);
		if(fatherlevel == 0) return TFather::SideNodeLocId(fatherside,node);
		if(fatherlevel == 1) return TFather::SideNodeLocId(fatherside,node-TFather::NSideNodes(fatherside))+TFather::NCornerNodes;//TFather::NSideNodes(fatherside);
		PZError << "Pr<TFather>::SideNodeLocId inconsistent side or node " << side
		<< ' ' << node << endl;
		return -1;
	}
	
	template<class TFather>
	int Pr<TFather>::SideDimension(int side) {
		if(side<0 || side >= NSides) {
			PZError << "Pr<TFather>::SideDimension side " << side << endl;
			return -1;
		}
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		if(ftns < 2) return TFather::SideDimension(fatherside);
		return TFather::SideDimension(fatherside)+1;
	}
	
	template<class TFather>
	void Pr<TFather>::HigherDimensionSides(int side, TPZStack<int> &high)
	{
		if(side <0 || side >= NSides) {
			PZError << "Pr<TFather>::HigherDimensionSides side "<< side << endl;
		}
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		int nhbefore = high.NElements();
		TFather::HigherDimensionSides(fatherside,high);
		int nhafter = high.NElements();
		int is;
		if(ftns < 2) 
		{
			high.Push(fatherside+2*TFather::NSides);
			for(is=nhbefore; is<nhafter; is++)
			{
				high.Push(high[is]+2*TFather::NSides);
				high[is] += ftns*TFather::NSides;
			}
		}
		if(ftns == 2)
		{
			for(is=nhbefore; is<nhafter; is++)
			{
				high[is] += ftns*TFather::NSides;
			}
		}
	}
	
	template<class TFather>
	void Pr<TFather>::CenterPoint(int side, TPZVec<REAL> &center) {
		//center.Resize(Dimension);
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		TFather::CenterPoint(fatherside,center);
		switch(ftns)
		{
			case 0:
				center[Dimension-1] = -1.;
				break;
			case 1:
				center[Dimension-1] = 1.;
				break;
			case 2:
			default:
				center[Dimension-1] = 0.;
		}
	}
	
	template<class TFather>
	TPZTransform Pr<TFather>::TransformElementToSide(int side){
		
		if(side<0 || side>=NSides)
		{
			PZError << "Pr<TFather>::TransformElementToSide called with side error\n";
			return TPZTransform(0,0);
		}
		
		int sidedim = SideDimension(side);
		TPZTransform t(sidedim,Dimension);//t(dimto,2)
		if(side == NSides) return t;
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		TPZTransform fathert = TFather::TransformElementToSide(fatherside);
		int fathersidedimension = sidedim;
		if(ftns == 2) fathersidedimension--;
		int id,jd;
		for(id=0; id<fathersidedimension; id++)
		{
			t.Sum()(id,0) = fathert.Sum()(id,0);
			for(jd=0; jd<TFather::Dimension; jd++)
			{
				t.Mult()(id,jd) = fathert.Mult()(id,jd);
			}
		}
		if(ftns == 2)
		{
			t.Mult()(sidedim-1,Dimension-1) = 1.;
		}
		return t;
	}
	
	template<class TFather>
	TPZTransform Pr<TFather>::TransformSideToElement(int side){
		
		if(side<0 || side>=NSides){
			PZError << "Pr<TFather>::TransformSideToElement side out range\n";
			return TPZTransform(0,0);
		}
		int sidedim = SideDimension(side);
		TPZTransform t(Dimension,sidedim);
		if(Dimension == sidedim) return t;
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		TPZTransform fathert = TFather::TransformSideToElement(fatherside);
		int fathersidedimension = sidedim;
		if(ftns == 2) fathersidedimension--;
		int id,jd;
		for(id=0; id<TFather::Dimension; id++)
		{
			t.Sum()(id,0) = fathert.Sum()(id,0);
			for(jd=0; jd<fathersidedimension; jd++)
			{
				t.Mult()(id,jd) = fathert.Mult()(id,jd);
			}
		}
		switch(ftns)
		{
			case 0:
				t.Sum()(Dimension-1,0) = -1.;
				break;
			case 1:
				t.Sum()(Dimension-1,0) = 1.;
				break;
			case 2:
				t.Mult()(Dimension-1,sidedim-1) = 1.;
				break;
		}
		return t;
		
	}
	
	template<class TFather>
	TPZTransform Pr<TFather>::SideToSideTransform(int sidefrom, int sideto)
	{
		if(sidefrom <0 || sidefrom >= NSides || sideto <0 || sideto >= NSides) {
			PZError << "Pr<TFather>::HigherDimensionSides sidefrom "<< sidefrom << 
			' ' << sideto << endl;
			return TPZTransform(0);
		}
		if(sidefrom == sideto) {
			return TPZTransform(SideDimension(sidefrom));
		}
		if(sidefrom == NSides-1) {
			return TransformElementToSide(sideto);
		}
		TPZStack<int> highsides;
		HigherDimensionSides(sidefrom,highsides);
		int nhigh = highsides.NElements();
		int is;
		for(is=0; is<nhigh; is++) {
			if(highsides[is] == sideto) {
				int dfr = SideDimension(sidefrom);
				int dto = SideDimension(sideto);
				TPZTransform trans(dto,dfr),toelement(Dimension,sidefrom),toside(sideto,Dimension);
				toelement = TransformSideToElement(sidefrom);
				toside = TransformElementToSide(sideto);
				trans = toside.Multiply(toelement);
				return trans;
			}
		}
		PZError << "Pr<TFather>::SideToSideTransform highside not found sidefrom "
		<< sidefrom << ' ' << sideto << endl;
		return TPZTransform(0);
	}
	
	template<class TFather>
	int Pr<TFather>::NumSides() {
		return TFather::NumSides()*3;
	}
	
	
	template<class TFather>
	int Pr<TFather>::NContainedSides(int side) {
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		if(ftns<2) return TFather::NContainedSides(fatherside);
		return TFather::NContainedSides(fatherside)*3;
	}
	
	template<class TFather>
	int Pr<TFather>::ContainedSideLocId(int side,int c) {
		int nsconnect = NContainedSides(side);
		if(c >= nsconnect)
		{
			PZError << "Pr<TFather>::ContainedSideLocId, side = " << side << " connect = " << c << endl;
			return -1;
		}
		int ftns = side/TFather::NSides;
		int fatherside = side%TFather::NSides;
		if(ftns == 0) return TFather::ContainedSideLocId(fatherside,c);
		if(ftns == 1) return TFather::ContainedSideLocId(fatherside,c)+TFather::NSides;
		int nsideconfather = TFather::NContainedSides(fatherside);
		int confather = c%nsideconfather;
		int level = c/nsideconfather;
		return TFather::ContainedSideLocId(fatherside,confather)+level*TFather::NSides;
	}
	
	template<class TFather>
	void Pr<TFather>::LowerDimensionSides(int side,TPZStack<int> &smallsides)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		int is;
		for(is=0; is<nsidecon-1; is++)
			smallsides.Push(ContainedSideLocId(side,is));
	}
	
	template<class TFather>
	void Pr<TFather>::LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget)
	{
		smallsides.Resize(0);
		int nsidecon = NContainedSides(side);
		for(int is = 0; is < nsidecon - 1; is++) {
			if (SideDimension(ContainedSideLocId(side,is)) == DimTarget) smallsides.Push(ContainedSideLocId(side,is));
		}
	}
	
	/**volume of the master element*/
	template<class TFather>
	REAL Pr<TFather>::RefElVolume()
	{
		return 2.0*TFather::RefElVolume();
	}
	
	template<class TFather>
	TPZIntPoints *Pr<TFather>::CreateSideIntegrationRule(int side, int order) {
		
		TPZIntPoints *integ = TFather::CreateSideIntegrationRule(side%TFather::NSides,order);
		if(side< 2*TFather::NSides) {
			return integ;
		}
		TPZIntPoints *integext = integ->PrismExtend(order);
		delete integ;
		return integext;
		
	}
	
	template<class TFather>
	std::string Pr<TFather>::StrType()
	{
		return "Pr:" + MElementType_Name(TFather::Type());
	}
	
	template<class TFather>
	std::string Pr<TFather>::StrType(int side)
	{
		if(side < 2*TFather::NSides)
		{
			return MElementType_Name(TFather::Type(side%(TFather::NSides)));
		}
		return "Pr:" + MElementType_Name(TFather::Type(side%(TFather::NSides)));
	}
	
	template<class TFather>
	bool Pr<TFather>::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide) {
		TPZTransform Transf = SideToSideTransform(NSides - 1, side);
		SidePar.Resize(SideDimension(side));
		Transf.Apply(InternalPar,SidePar);
		
		JacToSide = Transf.Mult();
		return true;
	}
	
	/**
	 * uses log4cxx to print the results of all methods
	 */
	template<class TFather>
	void Pr<TFather>::Diagnostic()
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << __PRETTY_FUNCTION__ ;
			sout << " Dimension : " << Dimension << " NCornerNodes " << NCornerNodes << 
			" NSides " << NSides;
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		{
			std::stringstream sout;
			sout << "Testing smallsides\n";
			TPZStack<int> smallsides;
			int is;
			for(is=0; is<NSides; is++)
			{
				LowerDimensionSides(is,smallsides);
				sout << " side " << is << " allsmallsides " << smallsides << endl;
				int id;
				for(id=0; id < SideDimension(is); id++)
				{
					smallsides.Resize(0);
					LowerDimensionSides(is,smallsides,id);
					sout << " side " << is << " targetdimension " << id << " smallsides " << smallsides << endl;
				}
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static void HigherDimensionSides(int side, TPZStack<int> &high);
		{
			std::stringstream sout;
			sout << "Testing HigherDimensionSides\n";
			int is ;
			for(is=0; is<NSides; is++)
			{
				TPZStack<int> high;
				HigherDimensionSides(is,high);
				sout << "Side " << is << " higherdimensionsides " << high << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static int NSideNodes(int side);
		{
			std::stringstream sout;
			sout <<  "Testing NSideNodes\n";
			int is; 
			for(is=0; is<NSides; is++)
			{
				sout << "Side " << is << " NSideNodes " << NSideNodes(is) << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		/**
		 * returns the local node number of the node "node" along side "side"
		 */
		//  static int SideNodeLocId(int side, int node);
		{
			std::stringstream sout;
			sout << "Testing SideNodeLocId\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				int nsnodes = NSideNodes(is);
				int sn;
				sout << "Side " << is << " sidenodelocid ";
				for(sn=0; sn<nsnodes; sn++)
				{
					sout << SideNodeLocId(is,sn) << " ";
				}
				sout << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static void CenterPoint(int side, TPZVec<REAL> &center);
		{
			std::stringstream sout;
			sout << "Testing CenterPoint\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				TPZManVector<REAL> center(Dimension);
				CenterPoint(is,center);
				sout << "Side " << is << center << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static REAL RefElVolume(){return 2.0*TFather::RefElVolume();}
		{
			std::stringstream sout;
			sout << "Testing RefelVolume\n";
			sout << "Volume " << RefElVolume() << endl;
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static int SideDimension(int side);
		{
			std::stringstream sout;
			sout << "Testing SideDimension\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				sout << "Side " << is << " sidedimension " << SideDimension(is) << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static TPZTransform SideToSideTransform(int sidefrom, int sideto);
		{
			std::stringstream sout;
			sout << "Testing SideToSideTransform\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				TPZStack<int> lower;
				LowerDimensionSides(is,lower);
				int nl = lower.NElements();
				int il;
				for(il=0; il<nl; il++)
				{
					TPZTransform tr;
					tr = SideToSideTransform(is,lower[il]);
					sout << "Transform from " << is << " to side " << lower[il] << " " << tr; 
				}
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static TPZTransform TransformElementToSide(int side);
		{
			std::stringstream sout;
			sout << "Testing TransformElementToSide\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				TPZTransform tr = TransformElementToSide(is);
				sout << "side " << is << " transform " << tr;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static TPZTransform TransformSideToElement(int side);
		{
			std::stringstream sout;
			sout << "Testing TransformSideToElement\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				sout << "side " << is << " transform " << TransformSideToElement(is);
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		{
			std::stringstream sout;
			sout << "Testing CreateSideIntegrationRule order 3\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				TPZIntPoints *integ = CreateSideIntegrationRule(is,3);
				sout << "side " << is << endl;
				integ->Print(sout);
				delete integ;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static std::string StrType() ;//{ return EOned;}
		{
			std::stringstream sout;
			sout << "Testing StrType\n";
			sout << StrType();
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static std::string StrType(int side);
		{
			std::stringstream sout;
			sout << "Testing StrType(side)\n";
			int is; 
			for(is=0; is<NSides; is++)
			{
				sout << "side " << is << " type " << StrType(is) << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static void MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);
		{
			std::stringstream sout;
			TPZManVector<REAL> par(Dimension,sqrt(2.));
			sout << "Testing MapToSide internalpar = " << par <<  "\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				TPZManVector<REAL> sidepar(SideDimension(is));
				TPZFMatrix<REAL> jactoside (SideDimension(is),Dimension);
				MapToSide(is,par,sidepar,jactoside);
				sout << "side " << is << " sidepar " << sidepar << " jactoside " << jactoside << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static int NSides;
		{
			std::stringstream sout;
			sout << "Testing NumSides\n";
			sout << "nsides " << NumSides();
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		//  static int NContainedSides(int side);
		{
			std::stringstream sout;
			sout << "Testing NContainedSides(side)\n";
			int is;
			for(is=0; is<NSides; is++)
			{
				sout << "side " << is << " nsideconnects " << NContainedSides(is) << endl;
			}
			LOGPZ_DEBUG(logger,sout.str());
		}
		//  static int ContainedSideLocId(int side, int c);
		{
			std::stringstream sout;
			sout << "Testing ContainedSideLocId(is,ic)\n";
			int is,ic;
			for(is=0; is<NSides; is++)
			{
				int nc = NContainedSides(is);
				sout << "side " << is << " connects ";
				for(ic=0; ic<nc; ic++)
				{
					sout << ContainedSideLocId(is,ic) << " ";
				}
				sout << endl;
			}
			
			LOGPZ_DEBUG(logger,sout.str());
		}
		
		
#endif
	}
	
	
	template<>
	REAL Pr< pztopology::TPZPoint >::RefElVolume()
	{
		return 2.0L;
	}
	
	template class Pr<TPZPoint>;
	template class Pr<Pr<TPZPoint> >;
	
}

