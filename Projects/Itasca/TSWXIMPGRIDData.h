//---------------------------------------------------------------------------

#ifndef TSWXIMPGRIDDataH
#define TSWXIMPGRIDDataH

#include <vector>
#include <map>
#include <string>
#include <iostream>

struct TSWXGridPoint{
  public:
  double fX[3];
  unsigned int fId;
  unsigned int fIndex;

  TSWXGridPoint()
  {
    //Do nothing
  }

  TSWXGridPoint(const TSWXGridPoint &cp){
    for(unsigned int i = 0; i < 3; i++) fX[i] = cp.fX[i];
    fId = cp.fId;
    fIndex = cp.fIndex;
  }

  TSWXGridPoint& operator=(const TSWXGridPoint &cp){
    for(unsigned int i = 0; i < 3; i++) fX[i] = cp.fX[i];
    fId = cp.fId;
    fIndex = cp.fIndex;
    return *this;
  }

  TSWXGridPoint(const unsigned int &index, const unsigned int &id, const double &x, const double &y, const double &z){
    fX[0] = x;
    fX[1] = y;
    fX[2] = z;
    fId = id;
    fIndex = index;
  }
};

template < int N >
class TSWXGridElement{
  public:
  unsigned int fIndices[N];
  unsigned int fId;
  unsigned int fIndex;

  TSWXGridElement(){
    fId = fIndex = -1;
    for(unsigned int i = 0; i < N; i++) fIndices[i] = -1;
  }

  TSWXGridElement(const TSWXGridElement<N> &cp){
    fIndex = cp.fIndex;
    fId = cp.fId;
    for(unsigned int i = 0; i < N; i++) fIndices[i] = cp.fIndices[i];
  }

  TSWXGridElement &operator=(const TSWXGridElement<N> &cp){
    fIndex = cp.fIndex;
    fId = cp.fId;
    for(unsigned int i = 0; i < N; i++) fIndices[i] = cp.fIndices[i];
    return *this;
  }

  TSWXGridElement(const unsigned int &index, const unsigned int &id, const std::vector<unsigned int> &indices){
    fIndex = index;
    fId = id;
    for(unsigned int i = 0; i < N; i++) fIndices[i] = indices[i];
  }
};

struct TSWXZoneGroup{
  public:
  std::string fGroupName;
  std::vector< unsigned int > fEls;

  TSWXZoneGroup()
  {
    //Do nothing
  }

  TSWXZoneGroup(const std::string &GroupName, const std::vector< unsigned int > &Els){
    fGroupName = GroupName;
    fEls = Els;
  }
  TSWXZoneGroup & operator=(const TSWXZoneGroup&cp){
    fGroupName = cp.fGroupName;
    fEls = cp.fEls;
    return *this;
  }
  TSWXZoneGroup(const TSWXZoneGroup&cp){
    fGroupName = cp.fGroupName;
    fEls = cp.fEls;
  }
};

class TSWXGridData{

  private:

//	CRITICAL_SECTION fGridCS, fGridElsCS, fZoneGrpCS;

  // fGridPointMapping[ id ] = index
  std::map< int, int > fGridPointMapping;

  // vector of gridpoints
  std::vector< TSWXGridPoint > fGridPoints;

  // fGridPointMapping[ id ] = index
  std::map< int, int > fB8ElementMapping, fW6ElementMapping;

  // vector of gridpoints
  std::vector< TSWXGridElement< 8 > > fB8Elements;

  std::vector< TSWXGridElement<6> > fW6Elements;

  std::vector< TSWXZoneGroup > fZoneGroups;

  public:

  TSWXGridData()
    : fGridPointMapping(), fGridPoints(),
      fB8ElementMapping(), fW6ElementMapping(), fB8Elements(),fW6Elements(){
//    InitializeCriticalSection(&fGridCS);
//    InitializeCriticalSection(&fGridElsCS);
//    InitializeCriticalSection(&fZoneGrpCS);
  }

  ~TSWXGridData(){
//    DeleteCriticalSection(&fGridCS);
//    DeleteCriticalSection(&fGridElsCS);
//    DeleteCriticalSection(&fZoneGrpCS);
  }

  unsigned int AddGridPoint( const unsigned int &id, const double &x, const double &y, const double &z){
//		EnterCriticalSection(&fGridCS);
    unsigned int index = this->fGridPoints.size();
    TSWXGridPoint pt(index,id,x,y,z);
    this->fGridPoints.push_back(pt);
    fGridPointMapping[ id ] = index;
//    LeaveCriticalSection(&fGridCS);
    return index;
  }

  unsigned int GetSizeGridPoint() const
  {
    return fGridPoints.size();
  }

  unsigned int AddB8Element( const unsigned int &id, const std::vector<unsigned int> & indices){
//		EnterCriticalSection(&fGridElsCS);
    unsigned int index = this->fB8Elements.size();
    TSWXGridElement<8> el(index,id,indices);
    this->fB8Elements.push_back(el);
    fB8ElementMapping[ id ] = index;
//    LeaveCriticalSection(&fGridElsCS);
    return index;
  }

  unsigned int AddW6Element( const unsigned int &id, const std::vector<unsigned int> & indices){
//		EnterCriticalSection(&fGridElsCS);
    unsigned int index = this->fW6Elements.size();
    TSWXGridElement<6> el(index,id,indices);
    this->fW6Elements.push_back(el);
    fW6ElementMapping[ id ] = index;
//    LeaveCriticalSection(&fGridElsCS);
    return index;
  }

  unsigned int GetSizeW6Element() const
  {
    return fW6Elements.size();
  }

  unsigned int GetSizeB8Element() const
  {
    return fB8Elements.size();
  }

  void AddZoneGroup(const std::string &GroupName, const std::vector< unsigned int > &Els){
//  	EnterCriticalSection(&fZoneGrpCS);
    TSWXZoneGroup zg(GroupName,Els);
    this->fZoneGroups.push_back(zg);
//    LeaveCriticalSection(&fZoneGrpCS);
  }

  unsigned int GetSizeZoneGroup() const
  {
    return fZoneGroups.size();
  }

  std::map< int, int > &GridPointMapping() { return fGridPointMapping; }

  const std::vector< TSWXGridPoint > &GridPoints() const { return fGridPoints; }

  std::map< int, int > &B8ElementMapping() { return fB8ElementMapping; }

  std::map< int, int > &W6ElementMapping() { return fW6ElementMapping; }

  const std::vector< TSWXGridElement<8> > &B8Elements() const{ return fB8Elements; }

  const std::vector< TSWXGridElement<6> > &W6Elements() const{ return fW6Elements; }

  const std::vector< TSWXZoneGroup > &ZoneGroups() const{ return fZoneGroups; }

  void Print(std::ostream &myfile);

  void Clear()
  {
    fGridPointMapping.clear();
    fGridPoints.resize(0);
    fB8ElementMapping.clear();
    fW6ElementMapping.clear();
    fB8Elements.resize(0);
    fW6Elements.resize(0);
    fZoneGroups.resize(0);
  }

};

//---------------------------------------------------------------------------
#endif
