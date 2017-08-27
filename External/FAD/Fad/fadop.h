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
#ifndef _fadop_h_
#define _fadop_h_

using namespace std;
//------------------------------- Fad binary operators ------------------------------------------


//------------------------------- Fad addition operators ------------------------------------------
template <typename L, typename R> class FadBinaryAdd {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryAdd() {}

  const L& left_; const R& right_;

public:
  FadBinaryAdd(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i) + right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i)+right_.fastAccessDx(i);}
};


template <typename L, typename R> class FadBinaryAdd<L, FadCst<R> > {
public:
  typedef typename L::value_type value_type_L;
  typedef R value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryAdd() {}

  const L& left_; const  FadCst<R> right_;

public:
  FadBinaryAdd(const L& left, const FadCst<R>& right) : left_(left), right_(right) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <typename L, typename R> class FadBinaryAdd< FadCst<L>, R> {
public:
  typedef L value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryAdd() {}

  const FadCst<L> left_; const R& right_;

public:
  FadBinaryAdd(const FadCst<L>& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  value_type dx(int i) const {return right_.dx(i);}
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i);}
};






//------------------------------- Fad substraction operators ------------------------------------------
template <typename L, typename R> class FadBinaryMinus {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryMinus() {}

  const L& left_; const R& right_;
  
public:
  FadBinaryMinus(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i) - right_.dx(i);}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) - right_.fastAccessDx(i);}
};


template <typename L, typename R> class FadBinaryMinus<L, FadCst<R> > {
public:
  typedef typename L::value_type value_type_L;
  typedef R value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryMinus() {}

  const L& left_; const FadCst<R> right_;

public:
  FadBinaryMinus(const L& left, const FadCst<R> & right) : left_(left), right_(right) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <typename L, typename R> class FadBinaryMinus< FadCst<L>, R> {
public:
  typedef L value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryMinus() {}

  const FadCst<L> left_; const R& right_;

public:
  FadBinaryMinus(const FadCst<L>& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return - right_.dx(i);}
  int size() const { return right_.size(); }

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - right_.fastAccessDx(i);}
};


//------------------------------- Fad multiplication operators ------------------------------------------
template <typename L, typename R> class FadBinaryMul {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
 protected:
  FadBinaryMul() {}

  const L& left_; const R& right_;

 public:
  FadBinaryMul(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val() + right_.dx(i) * left_.val();}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val() + right_.fastAccessDx(i) * left_.val();}

};

template <typename L, typename R> class FadBinaryMul<L, FadCst<R> > {
 public:
  typedef typename L::value_type value_type_L;
  typedef R value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
 protected:
  FadBinaryMul() {}

  const L& left_; const FadCst<R> right_;

 public:
  FadBinaryMul(const L& left, const FadCst<R>& right) : left_(left), right_(right) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val();}
};

template <typename L, typename R> class FadBinaryMul< FadCst<L>, R> {
 public:
  typedef L value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
 protected:
  FadBinaryMul() {}

  const FadCst<L> left_; const R& right_;

 public:
  FadBinaryMul(const FadCst<L>& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  right_.dx(i) * left_.val();}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i) * left_.val();}
};


//------------------------------- Fad division operators ------------------------------------------
template <typename L, typename R> class FadBinaryDiv {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:
  FadBinaryDiv() {}

  const L& left_; const R& right_;

 public:
  FadBinaryDiv(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (left_.dx(i) * right_.val() - right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (left_.fastAccessDx(i) * right_.val() - right_.fastAccessDx(i) * left_.val() ) 
					   / (right_.val() * right_.val()) ;}
};


template <typename L, typename R> class FadBinaryDiv<L, FadCst<R> > {
 public:
  typedef typename L::value_type value_type_L;
  typedef R value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:
  FadBinaryDiv() {}

  const L& left_; const FadCst<R> right_;

 public:
  FadBinaryDiv(const L& left, const FadCst<R>& right) : left_(left), right_(right) {;}
  ~FadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  left_.dx(i) / right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) / right_.val() ;}
};


template <typename L, typename R> class FadBinaryDiv< FadCst<L>, R> {
 public:
  typedef L value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:
  FadBinaryDiv() {}

  const FadCst<L> left_; const R& right_;

 public:
  FadBinaryDiv(const FadCst<L>& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryDiv() {;}

  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (- right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (- right_.fastAccessDx(i) * left_.val() ) 
					   / (right_.val() * right_.val()) ;}
};


//------------------------------- Fad pow function ------------------------------------------
template <typename L, typename R> class FadBinaryPow {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryPow() {}

  const L& left_; const R& right_;

public:
  FadBinaryPow(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return pow( left_.val(), right_.val() );}
  const value_type dx(int i) const 
    {
      return  (right_.dx(i)*log(left_.val())+right_.val()*left_.dx(i)/left_.val())
	*pow( left_.val(), right_.val() );
    }
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.fastAccessDx(i)*log(left_.val())+right_.val()*left_.fastAccessDx(i)/left_.val())
	*pow( left_.val(), right_.val() );
    }
};


template <typename L, typename R> class FadBinaryPow<L, FadCst<R> > {
public:
  typedef typename L::value_type value_type_L;
  typedef R value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryPow() {}

  const L& left_; const  FadCst<R> right_;

public:
  FadBinaryPow(const L& left, const FadCst<R> & right) : left_(left), right_(right) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return pow(left_.val(),right_.val()) ;}
  const value_type dx(int i) const 
    {
      return  (right_.val()*left_.dx(i)/left_.val())*pow( left_.val(), right_.val() );
    }
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.val()*left_.fastAccessDx(i)/left_.val())
	*pow( left_.val(), right_.val() );
    }
};


template <typename L, typename R> class FadBinaryPow< FadCst<L>, R> {
public:
  typedef L value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryPow() {}

  const FadCst<L> left_; const R& right_;

public:
  FadBinaryPow(const FadCst<L>& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryPow() {;}

  const value_type val() const {return pow(left_.val(),right_.val());}
  value_type dx(int i) const 
    {
      return (right_.dx(i)*log(left_.val()))*pow( left_.val(), right_.val() );
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.fastAccessDx(i)*log(left_.val()))
	*pow( left_.val(), right_.val() );
    }
};

template <typename L> class FadBinaryPow< L , int> {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<int> R;

protected:
  FadBinaryPow() {}

  const L& left_; const R right_;

public:
  FadBinaryPow(const L& left, const R& right) : left_(left), right_(right) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return pow(left_.val(),right_.val());}
  value_type dx(int i) const 
    {
      return right_.val()*pow( left_.val(), right_.val()-1);
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  right_.val() * pow( left_.val(), right_.val()-1 );
    }
};

#include "pzreal.h"
#include "fad.h"
//------------------------------- Fad operators ------------------------------------------
#define FAD_BIN_MACRO(OP, TYPE)                                                \
  /*A1 (a FadSuper object) vs A2 (another FadSuper object)*/                   \
  template <typename A1, typename A2,                                          \
            typename enable_if<is_convertible<A1*,FadSuper*>::value            \
                                 && is_convertible<A2*,FadSuper*>::value,      \
                                     int>::type * = nullptr>                   \
  inline FadExpr<TYPE<A1, A2>> OP(const A1 &v, const A2 &w) {    \
    typedef TYPE<A1, A2> expr_t;                             \
    return FadExpr<expr_t>(expr_t(v, w));                                      \
  }                                                                            \
  /*B1 (an arithmetic value) vs B2 (a FadSuper object) */                      \
  template <typename B1, typename B2,                                          \
            typename enable_if<(is_arithmetic_pz<B1>::value && \
                                     (is_convertible<B2*, FadSuper*>::value)), \
                                    int>::type * = nullptr>                    \
  inline FadExpr<TYPE<FadCst<B1>, B2>> OP(const B1 a, const B2 &e) {             \
    typedef TYPE<FadCst<B1>, B2> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(FadCst<B1>(a), e));                          \
  }                                                                            \
  /*C1 (a FadSuper object) vs C2 (an arithmetic value) */                                  \
  template <typename C1, typename C2,                                          \
            typename enable_if<((is_convertible<C1*, FadSuper*>::value) && \
                                     is_arithmetic_pz<C2>::value), \
                                    int>::type * = nullptr>                    \
  inline FadExpr<TYPE<C1, FadCst<C2>>> OP(const C1 &e, const C2 a) {                  \
    typedef TYPE<C1, FadCst<C2>> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(e, FadCst<C2>(a)));                          \
  }

FAD_BIN_MACRO(operator+,FadBinaryAdd)
FAD_BIN_MACRO(operator-,FadBinaryMinus)
FAD_BIN_MACRO(operator*,FadBinaryMul)
FAD_BIN_MACRO(operator/,FadBinaryDiv)

FAD_BIN_MACRO(pow,FadBinaryPow)

#undef FAD_BIN_MACRO

#endif
