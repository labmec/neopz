// Emacs will be in -*- Mode: c++ -*-
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//  Tiny benchmark of TinyFad with and without
//  exression template(ET)
// 
//********************************************************
//
//  This is a simple example to show how to use the
//  library
//
//********************************************************


#include <iostream>

// Fad with expression template
#include <Fad/fad.h>
// Fixed size Fad
#include <TinyFad/tinyfad.h>
// Fixed size Fad with expression template
#include <TinyFadET/tfad.h>

//*****************************************
//
// Compute the derivative of y = x*x with
// respect to x
//
//*****************************************
int main()
{
  

  float x = 2.f;
  const int nvar = 1;

   Fad<float> x1(x), y1;
  
  TinyFad<nvar,float> x2(x), y2;

  TFad<nvar,float> x3(x), y3;

  // activate xi variables :
  // the first argument is the order number
  // of the activated variable and the second
  // is the number of activated variables. It's
  // the same interface for the 3 classes
  x1.diff(0,nvar);

  x2.diff(0,nvar);

  x3.diff(0,nvar);

  TinyFad<2, double> variavel1, variavel2, y;

  variavel1 = 4.;
  variavel2 = 5.;



  variavel1.diff(0, 3);
  variavel2.diff(1,3);

  // funcao y = v1²*v2 + v2*3.;

  y = variavel1 * variavel1 * variavel2 + variavel2 * 3.;

  cout <<"teste" << y << endl << endl;

  // do computation
  y1 = x1*x1;

  y2 = x2*x2;

  y3 = x3*x3;

  cout << "Fad     : " << y1 << endl
       << "TinyFad : " << y2 << endl
       << "TFad    : " << y3 << endl;
  /*
  const int N =1;

  TFad<N,TFad<N,TFad<N,float> > > xthird, result;
  TFad<N,TFad<N,float> > xsecnd;
  TFad<N,float> xfirst(4.);

  xfirst.diff(0,N);

  xsecnd=xfirst;
  xsecnd.diff(0,N);

  xthird=xsecnd;
  xthird.diff(0,N);

  result = xthird*xthird*xthird;
  */

  const int N = 1;

  /*

  TFad<N,TFad<N,float> > xsecnd, result;
  TFad<N,float> xfirst(7.);

  xfirst.fastAccessDx(0)=5;

  xsecnd=xfirst;
  xsecnd.fastAccessDx(0)=5;

  result = xsecnd * xsecnd;

  cout << "\n\nTest: first(50x) and second (50) derivatives of (5x)^2 with x=.1.4\n" << result;
  cout << "\nx= " << xsecnd;


  result *= xsecnd;
  cout << "\n\nTest: first(3x^2) and second (6x) derivatives of x^3 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(4x^3) and second (12x^2x) derivatives of x^4 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(5x^4) and second (20x^3) derivatives of x^5 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(6x^5) and second (30x^4) derivatives of x^6 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;
  */
/*
  Fad<Fad<float> > xsecnd, result;
  Fad<float> xfirst(7.), init;
  Fad<Fad<Fad<float> > > divertido(1,0,Fad<Fad<float > > (1,0,Fad<float> (1,0,3))), resultado;

 cout << divertido;

  resultado = divertido * divertido;

  cout << "\ndivertido ao quadrado: " << resultado;


  resultado = resultado * divertido;

  cout << "\ndivertido ao cubo: " << resultado;

  resultado = resultado * divertido;

  cout << "\ndivertido ^4: " << resultado;*/
/*
  xfirst.diff(0, N);
  xfirst.fastAccessDx(0)=5.;

  init.diff(0,N);
  init = 5.;

  xsecnd=xfirst;
  xsecnd.diff(0,N);
  xsecnd.fastAccessDx(0)=init;

  result = xsecnd * xsecnd;

  cout << "\n\nTest: first(50x) and second (50) derivatives of (5x)^2 with x=.1.4\n" << result;
  cout << "\nx= " << xsecnd;


  result *= xsecnd;
  cout << "\n\nTest: first(3x^2) and second (6x) derivatives of x^3 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(4x^3) and second (12x^2x) derivatives of x^4 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(5x^4) and second (20x^3) derivatives of x^5 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;

  result *= xsecnd;
  cout << "\n\nTest: first(6x^5) and second (30x^4) derivatives of x^6 with x=4.\n" << result;
  cout << "\nx= " << xsecnd;
  */

  Fad<float> all(3,0,0);
  all.fastAccessDx(1)=1.;
  all.fastAccessDx(2)=1.;
  cout << "\n all differentiable = " << all;

 Fad<Fad<float> > teste(3,0,Fad<float> (3,0,7));


 Fad<Fad<float> > multipl(3,0, all), res;
 multipl.fastAccessDx(1)=1.;
 multipl.fastAccessDx(2)=1.;

 res *= 3.;

 res = teste * multipl;

 cout << "\nteste" << teste;

 cout << "\nmultipl" << multipl;

 cout << "\nteste * multipl" << res;

  return 1;
}
