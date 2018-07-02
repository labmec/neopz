#include "TSWXVector.h"
#include <fstream>

template class swx::vector<int>;

template class swx::vector<double>;

template class swx::vector< swx::vector<int> >;

template class swx::vector< swx::vector<double> >;

template class swx::vector< std::vector<double> >;

template class std::vector< swx::vector<double> >;
#if defined(SWX_BUILDER_2010) || defined (SWX_BUILDER_XE2)
template class swx::vector< System::UnicodeString >;

int teste(){
  swx::vector< System::UnicodeString > vec;
  std::ifstream file("c:\\Temp\\teste.txt");
  vec.Read(file);
  _di_IXMLNode myNode;
  vec.Write(myNode);
  vec.Read(myNode);
  return 1;
}
#endif



