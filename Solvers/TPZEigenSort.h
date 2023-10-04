#ifndef _TPZEIGENSORT_
#define _TPZEIGENSORT_
#include <pzerror.h>
#include <string>
#include <map>

//! Sorting method for calculated eigenvalues
enum class TPZEigenSort{
  Invalid,
  AbsAscending,/*!< Ascending magnitude*/
  AbsDescending,/*!< Descending magnitude*/
  RealAscending,/*!< Ascending real part*/
  RealDescending,/*!< Descending real part*/
  ImagAscending,/*!< Ascending imaginary part*/
  ImagDescending,/*!< Descending imaginary part*/
  TargetRealPart,/*!< Real part closest to target*/
  TargetImagPart,/*!< Imaginary part closest to target*/
  TargetMagnitude/*!< Magnitude closest to target*/
};

inline TPZEigenSort TPZEigenSortFromString(const std::string name){
  static const std::map<std::string, TPZEigenSort> stringvals {
    {"AbsAscending",TPZEigenSort::AbsAscending},
    {"AbsDescending",TPZEigenSort::AbsDescending},
    {"RealAscending",TPZEigenSort::RealAscending},
    {"RealDescending",TPZEigenSort::RealDescending},
    {"ImagAscending",TPZEigenSort::ImagAscending},
    {"ImagDescending",TPZEigenSort::ImagDescending},
    {"TargetRealPart",TPZEigenSort::TargetRealPart},
    {"TargetImagPart",TPZEigenSort::TargetImagPart},
    {"TargetMagnitude",TPZEigenSort::TargetMagnitude}
  };
  auto itr = stringvals.find(name);
  if( itr != stringvals.end() ) {
    return itr->second;
  }
  return TPZEigenSort::Invalid; 
  
}

inline std::string TPZEigenSortToString(const TPZEigenSort val){
  switch(val){
  case TPZEigenSort::Invalid: return "Invalid"; 
  case TPZEigenSort::AbsAscending: return "AbsAscending";
  case TPZEigenSort::AbsDescending: return "AbsDescending";
  case TPZEigenSort::RealAscending: return "RealAscending";
  case TPZEigenSort::RealDescending: return "RealDescending";
  case TPZEigenSort::ImagAscending: return "ImagAscending";
  case TPZEigenSort::ImagDescending: return "ImagDescending";
  case TPZEigenSort::TargetRealPart: return "TargetRealPart";
  case TPZEigenSort::TargetImagPart: return "TargetImagPart";
  case TPZEigenSort::TargetMagnitude: return "TargetMagnitude";
  }
}
#endif