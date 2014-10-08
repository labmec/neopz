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
//  KCC 3.3f (kai C++) and egcs 1.1.2 (cygnus)
//  are OK.
//  Memory Requirement : 
//   - 45 Mo with egcs 
//    (-Wall -O6 -mcpu=i686 -march=i686 -fstrict-aliasing)
//   - 32 Mo with KCC
//    (-w +K3 -O3)
//
//********************************************************

// ANSI C++ include
#include <iostream>
#include <iomanip>

// Fad includes
#include <TinyFadET/tfad.h>
#include <TinyFad/tinyfad.h>

// Include for time computation
#include <utils/timer.h>


using namespace std;

template <int Num> void Bench_4_Op(const TinyFad<Num> & x, const int nloop);


int main() {

  const int nloop = 1000000;//Increase or decrease 

  TinyFad<1> un; TinyFad<2> deux; TinyFad<3> trois;
  TinyFad<4> quatre; TinyFad<5> cinq; TinyFad<6> six;
  TinyFad<7> sept; TinyFad<8> huit; TinyFad<9> neuf;TinyFad<10> dix;
  TinyFad<11> onze; TinyFad<12> douze; TinyFad<13> treize;
  TinyFad<14> quatorze; TinyFad<15> quinze; TinyFad<16> seize;
  TinyFad<17> dixsept; TinyFad<18> dixhuit; TinyFad<19> dixneuf;TinyFad<20> vingt;


  cout << "Clean Cache" << endl;
  Bench_4_Op(neuf,nloop);

  cout << "_Size_____TFad______TinyFad_______Hand_______Member___\n" ;

  Bench_4_Op(un,nloop);
  Bench_4_Op(deux,nloop);
  Bench_4_Op(trois,nloop);
  Bench_4_Op(quatre,nloop);
  Bench_4_Op(cinq,nloop);
  Bench_4_Op(six,nloop);
  Bench_4_Op(sept,nloop);
  Bench_4_Op(huit,nloop);
  Bench_4_Op(neuf,nloop);
  Bench_4_Op(dix,nloop);
  Bench_4_Op(onze,nloop);
  Bench_4_Op(douze,nloop);
  Bench_4_Op(treize,nloop);
  Bench_4_Op(quatorze,nloop);
  Bench_4_Op(quinze,nloop);
  Bench_4_Op(seize,nloop);
  Bench_4_Op(dixsept,nloop);
  Bench_4_Op(dixhuit,nloop);
  Bench_4_Op(dixneuf,nloop);
  Bench_4_Op(vingt,nloop);
  
  return 0;
}





template <int Num> void Bench_4_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;
  //  int Num = x.size();

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  TFad<Num> x1(2.f), x2(3.f), x3(5.f), x4(0.f), x5(3.f);
  if (Num==1)
    x1.diff(0,Num);
  else {
    x1.diff(0,Num);
    x2.diff(1,Num);
  }

  for (j=0; j<ntest; j++){

    Timer mytime;
    mytime.start();

    for (i=0; i<nloop; ++i) 
      x4 = 2.f*(x1*x2)-x3+x3/x5+4.f;

    mytime.stop();
    time1 += mytime.elapsedSeconds();


   TinyFad<Num> y1(2.f), y2(3.f), y3(5.f), y4(0.f), y5(3.f);
   if (Num==1)
     y1.diff(0,Num);
   else {
     y1.diff(0,Num);
     y2.diff(1,Num);
   }


   mytime.start();

   for (i=0; i<nloop; ++i)
     y4 = 2.f*(y1*y2)-y3+y3/y5+4.f;

   mytime.stop();
   time2 += mytime.elapsedSeconds();


   float z1=2., z2=3., z3=5., z4, z5=3.;
   float *tz1, *tz2, *tz3, *tz4, *tz5;
   tz1 = new float[Num]; tz2 = new float[Num]; tz3 = new float[Num]; tz4 = new float[Num]; tz5 = new float[Num];
   for (k=0; k<Num; ++k) {
     tz5[k] = tz4[k] = tz3[k] = tz1[k] = tz2[k] = float(0.);
   }
   tz1[0]=1.;tz2[1]=1.;
   
   mytime.start();

   for (i=0; i<nloop; ++i) {
     z4 = 2.f*(z1*z2)-z3+z3/z5+4.f;
     for (k=0; k<Num; ++k)
       tz4[k] = 2.f*(z2*tz1[k]+z1*tz2[k])-tz3[k]+ (z5*tz3[k] - z3*tz5[k])/(z5*z5);
   }
   
   mytime.stop();
   time3 += mytime.elapsedSeconds();
   
   delete [] tz1; delete [] tz2; delete [] tz3; delete [] tz4; delete [] tz5;
   

   mytime.start();
   
   for (i=0; i<nloop; ++i) {
     x4.val() = 2.f*(x1.val()*x2.val())-x3.val()+x3.val()/x5.val()+4.f;
     for (k=0; k<Num; ++k)
       x4.fastAccessDx(k) =  2.f*(x1.fastAccessDx(k)*x2.val()+x1.val()*x2.fastAccessDx(k)) - x3.fastAccessDx(k) 
	 + (x5.val()*x3.fastAccessDx(k)-x3.val()*x5.fastAccessDx(k))/(x5.val()*x5.val());
   }
   
   mytime.stop();
   time4 += mytime.elapsedSeconds();


  }


  cout << setw(12) << time1 << setw(12) << time2 << setw(12) << time3 << setw(12) << time4 << endl;
}
