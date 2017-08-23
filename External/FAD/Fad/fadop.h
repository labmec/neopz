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

//------------------------------- Fad binary operators ------------------------------------------


//------------------------------- Fad addition operators ------------------------------------------
template <class L, class R> class FadBinaryAdd {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryAdd() {}

  const L& left_; const R& right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
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


template <class L> class FadBinaryAdd<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;
  
protected:
  FadBinaryAdd() {}

  const L& left_; const  R right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class FadBinaryAdd< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:
  FadBinaryAdd() {}

  const L left_; const R& right_;

public:
  FadBinaryAdd(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryAdd() {;}


  const value_type val() const {return left_.val() + right_.val();}
  value_type dx(int i) const {return right_.dx(i);}
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i);}
};






//------------------------------- Fad substraction operators ------------------------------------------
template <class L, class R> class FadBinaryMinus {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

protected:
  FadBinaryMinus() {}

  const L& left_; const R& right_;
  
public:
  FadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
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


template <class L> class FadBinaryMinus<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

protected:
  FadBinaryMinus() {}

  const L& left_; const R right_;

public:
  FadBinaryMinus(const L& left, const R & rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return left_.dx(i);}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i);}
};


template <class R> class FadBinaryMinus< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:
  FadBinaryMinus() {}

  const L left_; const R& right_;

public:
  FadBinaryMinus(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMinus() {;}


  const value_type val() const {return left_.val() - right_.val();}
  const value_type dx(int i) const {return - right_.dx(i);}
  int size() const { return right_.size(); }

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return - right_.fastAccessDx(i);}
};


//------------------------------- Fad multiplication operators ------------------------------------------
template <class L, class R> class FadBinaryMul {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
 protected:
  FadBinaryMul() {}

  const L& left_; const R& right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
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

template <class L> class FadBinaryMul<L, FadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;
  
 protected:
  FadBinaryMul() {}

  const L& left_; const R right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  left_.dx(i) * right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) * right_.val();}
};

template <class R> class FadBinaryMul< FadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;
  
 protected:
  FadBinaryMul() {}

  const L left_; const R& right_;

 public:
  FadBinaryMul(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryMul() {;}

  const value_type val() const {return left_.val() * right_.val() ;}
  const value_type dx(int i) const {return  right_.dx(i) * left_.val();}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return right_.fastAccessDx(i) * left_.val();}
};


//------------------------------- Fad division operators ------------------------------------------
template <class L, class R> class FadBinaryDiv {
 public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;
  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;

 protected:
  FadBinaryDiv() {}

  const L& left_; const R& right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
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


template <class L> class FadBinaryDiv<L, FadCst<typename L::value_type> > {
 public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;

 protected:
  FadBinaryDiv() {}

  const L& left_; const R right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryDiv() {;}


  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  left_.dx(i) / right_.val();}
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return left_.fastAccessDx(i) / right_.val() ;}
};


template <class R> class FadBinaryDiv< FadCst<typename R::value_type>, R> {
 public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

 protected:
  FadBinaryDiv() {}

  const L left_; const R& right_;

 public:
  FadBinaryDiv(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryDiv() {;}

  const value_type val() const {return left_.val() / right_.val();}
  const value_type dx(int i) const {return  (- right_.dx(i) * left_.val() ) / (right_.val() * right_.val()) ;}
  int size() const { return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const { return (- right_.fastAccessDx(i) * left_.val() ) 
					   / (right_.val() * right_.val()) ;}
};


//------------------------------- Fad pow function ------------------------------------------
template <class L, class R> class FadBinaryPow {
public:
  typedef typename L::value_type value_type_L;
  typedef typename R::value_type value_type_R;

  typedef typename NumericalTraits<value_type_L,value_type_R>::promote value_type;
  
protected:
  FadBinaryPow() {}

  const L& left_; const R& right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow( left_.val(), right_.val() );}
  const value_type dx(int i) const 
    {
      return  (right_.dx(i)*std::log(left_.val())+right_.val()*left_.dx(i)/left_.val())
	*std::pow( left_.val(), right_.val() );
    }
  int size() const {
    int lsz = left_.size(), rsz = right_.size();
    return max(lsz, rsz);
  }

  bool hasFastAccess() const { return left_.hasFastAccess() && right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.fastAccessDx(i)*std::log(left_.val())+right_.val()*left_.fastAccessDx(i)/left_.val())
	*std::pow( left_.val(), right_.val() );
    }
};


template <class L> class FadBinaryPow<L, FadCst<typename L::value_type> > {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<value_type> R;
  
protected:
  FadBinaryPow() {}

  const L& left_; const  R right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow(left_.val(),right_.val()) ;}
  const value_type dx(int i) const 
    {
      return  (right_.val()*left_.dx(i)/left_.val())*std::pow( left_.val(), right_.val() );
    }
  int size() const { return left_.size();}

  bool hasFastAccess() const { return left_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.val()*left_.fastAccessDx(i)/left_.val())
	*std::pow( left_.val(), right_.val() );
    }
};


template <class R> class FadBinaryPow< FadCst<typename R::value_type>, R> {
public:
  typedef typename R::value_type value_type;
  typedef FadCst<value_type> L;

protected:
  FadBinaryPow() {}

  const L left_; const R& right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}

  const value_type val() const {return std::pow(left_.val(),right_.val());}
  value_type dx(int i) const 
    {
      return (right_.dx(i)*std::log(left_.val()))*std::pow( left_.val(), right_.val() );
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  (right_.fastAccessDx(i)*std::log(left_.val()))
	*std::pow( left_.val(), right_.val() );
    }
};

template <class L> class FadBinaryPow< L , int> {
public:
  typedef typename L::value_type value_type;
  typedef FadCst<int> R;

protected:
  FadBinaryPow() {}

  const L& left_; const R right_;

public:
  FadBinaryPow(const L& left, const R& rigth) : left_(left), right_(rigth) {;}
  ~FadBinaryPow() {;}


  const value_type val() const {return std::pow(left_.val(),right_.val());}
  value_type dx(int i) const 
    {
      return right_.val()*std::pow( left_.val(), right_.val()-1);
    }
  int size() const {return right_.size();}

  bool hasFastAccess() const { return right_.hasFastAccess();}
  value_type fastAccessDx(int i) const 
    {
      return  right_.val() * std::pow( left_.val(), right_.val()-1 );
    }
};

#include "pzreal.h"
//------------------------------- Fad operators ------------------------------------------
#define FAD_BIN_MACRO(OP, TYPE)                                                \
  /*FadExpr(A1) vs FadExpr(A2)*/                                               \
  template <typename A1, typename A2>                                          \
  inline FadExpr<TYPE<FadExpr<A1>, FadExpr<A2>>> OP(const FadExpr<A1> &v,      \
                                                    const FadExpr<A2> &w) {    \
    typedef TYPE<FadExpr<A1>, FadExpr<A2>> expr_t;                             \
    return FadExpr<expr_t>(expr_t(v, w));                                      \
  }                                                                            \
  /*FadExpr(B1) vs FadCst(B1)*/                                                \
  template <typename B1>                                                       \
  inline FadExpr<TYPE<FadExpr<B1>, FadCst<typename B1::value_type>>> OP(       \
      const FadExpr<B1> &e, const typename B1::value_type &t) {                \
    typedef typename B1::value_type A;                                         \
    typedef TYPE<FadExpr<B1>, FadCst<A>> expr_t;                               \
    return FadExpr<expr_t>(expr_t(e, FadCst<A>(t)));                           \
  }                                                                            \
  /*FadCst(C1) vs Fad(C1)*/                                                    \
  template <typename C1>                                                       \
  inline FadExpr<TYPE<FadCst<C1>, Fad<C1>>> OP(const C1 &a,                    \
                                               const Fad<C1> &e) {             \
    typedef TYPE<FadCst<C1>, Fad<C1>> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(FadCst<C1>(a), e));                          \
  }                                                                            \
  /*Fad(D1) vs FadCst(D1)*/                                                    \
  template <typename D1>                                                       \
  inline FadExpr<TYPE<Fad<D1>, FadCst<D1>>> OP(const Fad<D1> &e,               \
                                               const D1 &a) {                  \
    typedef TYPE<Fad<D1>, FadCst<D1>> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(e, FadCst<D1>(a)));                          \
  }                                                                            \
  /*FadCst(E1) vs Fad(E2), E1 is arithmetic*/                                  \
  template <typename E1, typename E2,                                          \
            typename std::enable_if<(std::is_integral<E1>::value ||            \
                                     is_complex_or_floating_point<E1>::value), \
                                    int>::type * = nullptr>                    \
  inline FadExpr<TYPE<FadCst<E1>, Fad<E2>>> OP(const E1 &a,                    \
                                               const Fad<E2> &e) {             \
    typedef TYPE<FadCst<E1>, Fad<E2>> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(FadCst<E1>(a), e));                          \
  }                                                                            \
  /*Fad(F1) vs FadCst(F2), F2 is arithmetic*/                                  \
  template <typename F1, typename F2,                                          \
            typename std::enable_if<(std::is_integral<F2>::value ||            \
                                     is_complex_or_floating_point<F2>::value), \
                                    int>::type * = nullptr>                    \
  inline FadExpr<TYPE<Fad<F1>, FadCst<F2>>> OP(const Fad<F1> &e,               \
                                               const F2 &a) {                  \
    typedef TYPE<Fad<F1>, FadCst<F2>> expr_t;                                  \
    return FadExpr<expr_t>(expr_t(e, FadCst<F2>(a)));                          \
  }                                                                            \
  /*FadCst(G1::value_type) vs FadExpr(G1)*/                                    \
  template <typename G1>                                                       \
  inline FadExpr<TYPE<FadCst<typename G1::value_type>, FadExpr<G1>>> OP(       \
      const typename G1::value_type &t, const FadExpr<G1> &e) {                \
    typedef typename G1::value_type A;                                         \
    typedef TYPE<FadCst<A>, FadExpr<G1>> expr_t;                               \
    return FadExpr<expr_t>(expr_t(FadCst<A>(t), e));                           \
  }                                                                            \
  /*FadExpr(H1) vs Fad(H1::value_type)*/                                       \
  template <typename H1>                                                       \
  inline FadExpr<TYPE<FadExpr<H1>, Fad<typename H1::value_type>>> OP(          \
      const FadExpr<H1> &e, const Fad<typename H1::value_type> &v) {           \
    typedef TYPE<FadExpr<H1>, Fad<typename H1::value_type>> expr_t;            \
    return FadExpr<expr_t>(expr_t(e, v));                                      \
  }                                                                            \
  /*Fad(I1) vs Fad(I1)*/                                                       \
  template <typename I1>                                                       \
  inline FadExpr<TYPE<Fad<I1>, Fad<I1>>> OP(const Fad<I1> &e1,                 \
                                            const Fad<I1> &e2) {               \
    typedef TYPE<Fad<I1>, Fad<I1>> expr_t;                                     \
    return FadExpr<expr_t>(expr_t(e1, e2));                                    \
  }                                                                            \
  /*Fad(J1::value_type) vs FadExpr(J1)*/                                       \
  template <typename J1>                                                       \
  inline FadExpr<TYPE<Fad<typename J1::value_type>, FadExpr<J1>>> OP(          \
      const Fad<typename J1::value_type> &v, const FadExpr<J1> &e) {           \
    typedef TYPE<Fad<typename J1::value_type>, FadExpr<J1>> expr_t;            \
    return FadExpr<expr_t>(expr_t(v, e));                                      \
  }




FAD_BIN_MACRO(operator+,FadBinaryAdd)
FAD_BIN_MACRO(operator-,FadBinaryMinus)
FAD_BIN_MACRO(operator*,FadBinaryMul)
FAD_BIN_MACRO(operator/,FadBinaryDiv)

FAD_BIN_MACRO(pow,FadBinaryPow)

#undef FAD_BIN_MACRO


#endif
