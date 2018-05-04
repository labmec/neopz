// Emacs will be in -*- Mode: c++ -*-
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses, 
//         templates : new C++ techniques 
//            for scientific computing 
// 
//********************************************************
//
//  A short implementation ( not all operators and 
//  functions are overloaded ) of 1st order Automatic
//  Differentiation in forward mode (FAD) using
//  EXPRESSION TEMPLATES.
//
//********************************************************
#ifndef _fad_h_
#define _fad_h_

#include <algorithm>

#include <cmath>
#include <math.h>

#include <utils/vectors.h>
#include <utils/promote.h>
#include <type_traits>

#include <Hash/TPZHash.h>

using namespace std;

template <class T> class FadExpr;
template <class T> class FadCst;

template <class T> class FadUnaryPlus;
template <class T> class FadUnaryMin;

template <class L, class R> class FadBinaryAdd;
template <class L, class R> class FadBinaryMinus;
template <class L, class R> class FadBinaryMul;
template <class L, class R> class FadBinaryDiv;

class FadSuper {
protected :
    FadSuper(){}
    FadSuper(const FadSuper&){};
    ~FadSuper(){};
};

template <class T> class Fad : public FadSuper {
public:
  typedef T value_type;
  
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  inline void copy(const Fad<R>& rhs)
  {
    const Vector<R> &xdx = rhs.dx();
    dx_.resize(xdx.size());
    size_t sz = rhs.dx().size();
    for (size_t i = 0; i < sz; ++i) {
        dx_[i] = xdx[i];
    }
    val_ = rhs.val();
  }
  
protected:
  T val_;
  Vector<T> dx_;

public:

  T defaultVal;

  Fad() : val_( T(0.f)), dx_(), defaultVal(T(0)) {;}
  
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  Fad(const R & x) : val_(x), dx_(), defaultVal(T(0)) {;}
  Fad(const int sz, const T & x) : val_(x), dx_(sz,T(0)), defaultVal(T(0)) {;}
  Fad(const int sz, const int i, const T & x) : val_(x), dx_(sz,T(0)), defaultVal(T(0))
    {dx_[i]=1.;}
  Fad(const int sz, const T & x, const T & dx) : val_(x), dx_(sz, dx), defaultVal(T(0)) {;}
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  Fad(const Fad<R> & rhs) : defaultVal(T(0))
  {
    copy(rhs);
  }
  //template <class ExprT, typename enable_if<is_convertible<typename ExprT::value_type,T>::value,int>::type * = nullptr> 
  template <class ExprT> 
    inline Fad(const FadExpr<ExprT>& fadexpr) : 
    val_(fadexpr.val()), 
    dx_(fadexpr.size()), 
    defaultVal(T(0))
    {
      int sz = fadexpr.size();

      if ( sz ) {
        for(int i=0; i<sz; ++i) 
          dx_[i] = fadexpr.dx(i);
      }
    }

  ~Fad(){;}

  void diff(const int ith, const int n);

  const Vector<T>& dx() const { return dx_;}

  const T& val()     const { return val_;}
  T& val()    { return val_;}

  bool hasFastAccess() const { return dx_.size()!=0;}

  T& fastAccessDx(int i) { return dx_[i];}
  const T& fastAccessDx(int i) const { return dx_[i];}
  const T& d(int i) const { return dx_[i];}
/*
  const T dx(int i) const { T tmp= dx_.size()? dx_[i]:T(0); return tmp;}
  */

  const T& dx(int i) const { return dx_.size() ? dx_[i] : defaultVal;}
  T dx(int i) { T tmp= dx_.size()? dx_[i]:defaultVal; return tmp;}

  int size() const { return dx_.size();}

  template <class ExprT> Fad<T> & operator=(const FadExpr<ExprT>& fadexpr);

  FadExpr< FadUnaryPlus< Fad<T> > > operator+ () const;
  FadExpr< FadUnaryMin< Fad<T> > > operator- () const;

  Fad<T>& operator+= (const T& x);
  Fad<T>& operator-= (const T& x);
  Fad<T>& operator*= (const T& x);
  Fad<T>& operator/= (const T& x);

  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  inline Fad<T>& operator+= (const Fad<R>& x){
    int xsz = x.size(), sz = dx_.size();

    if (xsz) {
      R* xdx = x.dx().begin();
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] += xdx[i];
      }
      else {
        dx_.resize(xsz);
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<xsz; ++i)
          dxp[i] = xdx[i];
      }
    }

    val_ += x.val();

    return *this;
  }
  
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  inline Fad<T>& operator-= (const Fad<R>& x){
    int xsz = x.size(), sz = dx_.size();

    if (xsz) {
      R* RESTRICT xdx = x.dx().begin();
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] -= xdx[i];
      }
      else {
        dx_.resize(xsz);
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<xsz; ++i)
          dxp[i] = - xdx[i];
      }
    }

    val_ -= x.val();

    return *this;
  }
  
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  inline Fad<T>& operator*= (const Fad<R>& x){
    int xsz = x.size(), sz = dx_.size();
    R xval = x.val();

    if (xsz) {
      R* RESTRICT xdx = x.dx().begin();
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] = val_ * xdx[i] + dxp[i] * xval;
      }
      else {
        dx_.resize(xsz);
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<xsz; ++i)
          dxp[i] = val_ * xdx[i];
      }
    }
    else {
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] *= xval;
      }
    }

    val_ *= xval;

    return *this;
  }
  
  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
  inline Fad<T>& operator/= (const Fad<R>& x){
    int xsz = x.size(), sz = dx_.size();
    R xval = x.val();

    if (xsz) {
      R* RESTRICT xdx = x.dx().begin(); 
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] = (dxp[i]*xval - val_*xdx[i])/ (xval*xval);
      }
      else {
        dx_.resize(xsz);
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<xsz; ++i)
          dxp[i] = -val_ * xdx[i] / (xval*xval);
      }
    }
    else {
      if (sz) {
        T* RESTRICT dxp = dx_.begin();
        for (int i=0; i<sz; ++i)
          dxp[i] /= xval;
      }
    }

    val_ /= x.val();

    return *this;
  }

  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
    inline Fad<T> &operator=(const R& val) 
    {
      val_ = val;

      if ( dx_.size() ) dx_ = T(0.);

      return *this;
    }

  template<typename R, typename enable_if<is_convertible<R,T>::value,int>::type * = nullptr>
    inline Fad<T> &operator=(const Fad<R>& rhs) 
    {
      if ( this != (Fad<T>*)(&rhs) ) copy(rhs);

      return *this;
    }
  
  virtual int ClassId() const;
  
  template <class ExprT> Fad<T>& operator*= (const FadExpr<ExprT>& fadexpr);
  template <class ExprT> Fad<T>& operator/= (const FadExpr<ExprT>& fadexpr);
  template <class ExprT> Fad<T>& operator+= (const FadExpr<ExprT>& fadexpr);
  template <class ExprT> Fad<T>& operator-= (const FadExpr<ExprT>& fadexpr);
    
  friend ostream& operator<< (ostream& stream, const Fad<T>& x)
  {
      stream << x.val();
//      for (int i=0; i<x.dx_.size(); i++) {
//          stream << x.dx_[i] << " ";
//      }
      return stream;
  }
    
  friend istream& operator>> (istream& stream, Fad<T>& x)
  {
      return stream >> x.val();
  }
    
};

template<class T> inline
T fabs(const Fad<T> &val)
{
    return fabs(val.val());
}

//template<class T> inline
//double Norm(const Fad<T> &val)
//{
//    DebugStop();
//}


template <class T> inline  void Fad<T>::diff(const int ith, const int n) 
{ 
  if ( dx_.empty() ) dx_.resize(n);

  dx_ = T(0.);
  dx_[ith] = T(1.);

}

template <class T> template <class ExprT> inline Fad<T> & Fad<T>::operator=(const FadExpr<ExprT>& fadexpr) 
{
  int sz = fadexpr.size();

  if ( sz != dx_.size() ) dx_.resize(sz);

  if ( sz ) {
    T* RESTRICT dxp = dx_.begin();
    if (fadexpr.hasFastAccess())
      for(int i=0; i<sz; ++i)
	dxp[i] = fadexpr.fastAccessDx(i);
    else
      for(int i=0; i<sz; ++i)
	dxp[i] = fadexpr.dx(i);
  }
  
  val_ = fadexpr.val();
  
  return *this;
}

template <class T> inline FadExpr< FadUnaryPlus< Fad<T> > >
Fad<T>::operator+ () const
{
  return FadExpr< FadUnaryPlus< Fad<T> > >(*this);
}

template <class T> inline FadExpr< FadUnaryMin< Fad<T> > >
Fad<T>::operator- () const
{
  return FadExpr< FadUnaryMin< Fad<T> > >(*this);
}


template <class T> inline  Fad<T> & Fad<T>::operator+= (const T& val)
{
  val_ += val;

  return *this;
}

template <class T> inline  Fad<T> & Fad<T>::operator-= (const T& val)
{
  val_ -= val;

  return *this;
}

template <class T> inline  Fad<T> & Fad<T>::operator*= (const T& val)
{
  val_ *= val;

  int sz = dx_.size();
  if ( sz ) {
    T* RESTRICT dxp = dx_.begin();
    for (int i=0; i<sz;++i)
      dxp[i] *= val;
  }

  return *this;
}

template <class T> inline  Fad<T> & Fad<T>::operator/= (const T& val)
{
  val_ /= val;

  int sz = dx_.size();
  if ( sz ) {
    T* RESTRICT dxp = dx_.begin();
    for (int i=0; i<sz;++i)
      dxp[i] /= val;
  }

  return *this;
}



template <class T> template <class ExprT> inline  Fad<T> & Fad<T>::operator+= (const FadExpr<ExprT>& x)
{
  int xsz = x.size(), sz = dx_.size();

  if (xsz) {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for (int i=0; i<sz; ++i)
	  dxp[i] += x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  dxp[i] += x.dx(i);
    }
    else {
      dx_.resize(xsz);
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for (int i=0; i<xsz; ++i)
	  dxp[i] = x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  dxp[i] = x.dx(i);
    }
  }

  val_ += x.val();

  return *this;
}

template <class T> template <class ExprT> inline  Fad<T> & Fad<T>::operator-= (const FadExpr<ExprT>& x)
{
  int xsz = x.size(), sz = dx_.size();

  if (xsz) {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  dxp[i] -= x.fastAccessDx(i);
      else
	for (int i=0; i<sz; ++i)
	  dxp[i] -= x.dx(i);
    }
    else {
      dx_.resize(xsz);
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  dxp[i] = -x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  dxp[i] = -x.dx(i);
    }
  }

  val_ -= x.val();


  return *this;
}

template <class T> template <class ExprT> inline  Fad<T> & Fad<T>::operator*= (const FadExpr<ExprT>& x)
{
  int xsz = x.size(), sz = dx_.size();
  T xval = x.val();

  if (xsz) {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  dxp[i] = val_ * x.fastAccessDx(i) + dxp[i] * xval;
      else
	for (int i=0; i<sz; ++i)
	  dxp[i] = val_ * x.dx(i) + dxp[i] * xval;
    }
    else {
      dx_.resize(xsz);
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  dxp[i] = val_ * x.fastAccessDx(i);
      else
	for (int i=0; i<xsz; ++i)
	  dxp[i] = val_ * x.dx(i);
    }
  }
  else {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      for (int i=0; i<sz; ++i)
	dxp[i] *= xval;
    }
  }

  val_ *= xval;

  return *this;
}

template <class T> template <class ExprT> inline  Fad<T> & Fad<T>::operator/= (const FadExpr<ExprT>& x)
{
  int xsz = x.size(), sz = dx_.size();
  T xval = x.val();

  if (xsz) {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<sz; ++i)
	  dxp[i] = ( dxp[i]*xval - val_*x.fastAccessDx(i) )/ (xval*xval);
      else
	for (int i=0; i<sz; ++i)
	  dxp[i] = ( dxp[i]*xval - val_*x.dx(i) )/ (xval*xval);
    }
    else {
      dx_.resize(xsz);
      T* RESTRICT dxp = dx_.begin();
      if (x.hasFastAccess())
	for(int i=0; i<xsz; ++i)
	  dxp[i] = - val_*x.fastAccessDx(i) / (xval*xval);
      else
	for (int i=0; i<xsz; ++i)
	  dxp[i] = -val_ * x.dx(i) / (xval*xval);
    }
  }
  else {
    if (sz) {
      T* RESTRICT dxp = dx_.begin();
      for (int i=0; i<sz; ++i)
	dxp[i] /= xval;
    }
  }

  val_ /= xval;

  return *this;
}




//------------------------------- Fad ostream operator ------------------------------------------
template <class T> inline ostream& operator << (ostream& os, const Fad<T>& a)
{
  os.setf(ios::fixed,ios::floatfield);
  os.width(12);
  os << a.val() << "  [";

  
  for (int i=0; i< a.dx().size(); i++) {
     os.width(12);
     os << a.dx(i);
  }

  os << "]\n";
  return os;
}

//------------------------------- Fad expression ------------------------------------------
template < class T > class FadExpr : public FadSuper {
public:
  typedef typename T::value_type value_type;

protected:
  FadExpr() {}

  T fadexpr_;

public:
  explicit FadExpr(const T& fadexpr) : fadexpr_(fadexpr) {;}

  value_type val()     const { return fadexpr_.val();}
  value_type dx(int i) const { return fadexpr_.dx(i);}
  int size() const {return fadexpr_.size();}

  bool hasFastAccess() const { return fadexpr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return fadexpr_.fastAccessDx(i);}
};

//------------------------------- Fad constant ------------------------------------------
template <class T> class FadCst : public FadSuper {
public:
    typedef T value_type;
    
protected:
    FadCst() {}
    
    const T constant_;
    const T defaultVal;
    
public:
    explicit FadCst(const T &value) : constant_(value), defaultVal(T(0)) {}
    
    const value_type &val() const { return constant_; }
    const value_type dx(int i) const { return defaultVal; }
    int size() const { return 0; }
    
    bool hasFastAccess() const { return 1; }
    const value_type &fastAccessDx(int i) const {
        return defaultVal;
    } // any i-th derivative of a constant is zero (i>0)
};

//------------------------------- Fad unary + ------------------------------------------
template < class T > class FadUnaryPlus {
public:
  typedef typename T::value_type value_type;

protected:
  FadUnaryPlus() {}

  const T& expr_;

public:
  FadUnaryPlus(const T& value) : expr_(value) {;}

  const value_type val()     const { return expr_.val();}
  const value_type dx(int i) const { return expr_.dx(i);}
  int size() const {return expr_.size();}
  
  bool hasFastAccess() const { return expr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return expr_.fastAccessDx(i);}
};

//------------------------------- Fad unary - ------------------------------------------
template < class T > class FadUnaryMin {
public:
  typedef typename T::value_type value_type;

protected:
  FadUnaryMin() {}

  const T& expr_;

public:
  FadUnaryMin(const T& value) : expr_(value) {;}

  const value_type val()     const { return - expr_.val();}
  const value_type dx(int i) const { return - expr_.dx(i);}
  int size() const {return expr_.size();}
  
  bool hasFastAccess() const { return expr_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - expr_.fastAccessDx(i);}
};

template <class T> inline
FadExpr< FadUnaryPlus< FadExpr<T> > >
operator + (const FadExpr<T>& expr)
{
  typedef FadUnaryPlus< FadExpr<T> > expr_t;

  return FadExpr< expr_t >( expr_t(expr) );
}

template <class T> inline
FadExpr< FadUnaryMin< FadExpr<T> > >
operator - (const FadExpr<T>& expr)
{
  typedef FadUnaryMin< FadExpr<T> > expr_t;

  return FadExpr< expr_t >( expr_t(expr) );
}

template <class T> 
int Fad<T>::ClassId() const{
    return Hash("Fad") ^ (ClassIdOrHash<T>()<<1);
}


#include <Fad/fadlog.h>
#include <Fad/fadop.h>
#include <Fad/fadfunc.h>

#endif
