// Emacs will be in -*- Mode: c++ -*-
//
// ********** DO NOT REMOVE THIS BANNER **********
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//  Big benchmark of TinyFad with and without
//  exression template(ET)
// 
//********************************************************
//
//  KCC 3.3f (kai C++) and egcs 1.1.2 (cygnus)
//  are OK.
//  Memory Requirement : 
//   - 108 Mo with egcs 
//    (-Wall -O6 -mcpu=i686 -march=i686 -fstrict-aliasing)
//   - 83 Mo with KCC
//    (-w +K3 -O3)
//
//********************************************************

// ANSI C++ include
#include <iostream>
#include <iomanip>

// Include for time computation
#include <utils/timer.h>

// Fad includes
#include <TinyFadET/tfad.h>
#include <TinyFad/tinyfad.h>

using namespace std;

template <int Num> void Bench_1_Op(const TinyFad<Num> & x, const int nloop);

template <int Num> void Bench_2_Op(const TinyFad<Num> & x, const int nloop);

template <int Num> void Bench_3_Op(const TinyFad<Num> & x, const int nloop);

template <int Num> void Bench_4_Op(const TinyFad<Num> & x, const int nloop);

template <int Num> void Bench_10_Op(const TinyFad<Num> & x, const int nloop);

template <int Num> void Bench_20_Op(const TinyFad<Num> & x, const int nloop);


int main() {

  int nloop = 1000000;//100
  TinyFad<1> un; TinyFad<2> deux; TinyFad<3> trois;
  TinyFad<4> quatre; TinyFad<5> cinq; TinyFad<6> six;
  TinyFad<7> sept; TinyFad<8> huit; TinyFad<9> neuf;TinyFad<10> dix;
  TinyFad<11> onze; TinyFad<12> douze; TinyFad<13> treize;
  TinyFad<14> quatorze; TinyFad<15> quinze; TinyFad<16> seize;
  TinyFad<17> dixsept; TinyFad<18> dixhuit; TinyFad<19> dixneuf;TinyFad<20> vingt;

  cout << "Clean Cache" << endl;
  Bench_4_Op(neuf,nloop);

  cout << "_Size_____TFad______TinyFad_______Hand_______Member___\n" ;

//   Bench_1_Op(un,nloop);
//   Bench_1_Op(deux,nloop);
//   Bench_1_Op(trois,nloop);
//   Bench_1_Op(quatre,nloop);
//   Bench_1_Op(cinq,nloop);
//   Bench_1_Op(six,nloop);
//   Bench_1_Op(sept,nloop);
//   Bench_1_Op(huit,nloop);
//   Bench_1_Op(neuf,nloop);
//   Bench_1_Op(dix,nloop);

//   Bench_2_Op(un,nloop);
//   Bench_2_Op(deux,nloop);
//   Bench_2_Op(trois,nloop);
//   Bench_2_Op(quatre,nloop);
//   Bench_2_Op(cinq,nloop);
//   Bench_2_Op(six,nloop);
//   Bench_2_Op(sept,nloop);
//   Bench_2_Op(huit,nloop);
//   Bench_2_Op(neuf,nloop);
//   Bench_2_Op(dix,nloop);

//   Bench_3_Op(un,nloop);
//   Bench_3_Op(deux,nloop);
//   Bench_3_Op(trois,nloop);
//   Bench_3_Op(quatre,nloop);
//   Bench_3_Op(cinq,nloop);
//   Bench_3_Op(six,nloop);
//   Bench_3_Op(sept,nloop);
//   Bench_3_Op(huit,nloop);
//   Bench_3_Op(neuf,nloop);
//   Bench_3_Op(dix,nloop);

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

  Bench_10_Op(un,nloop);
  Bench_10_Op(deux,nloop);
  Bench_10_Op(trois,nloop);
  Bench_10_Op(quatre,nloop);
  Bench_10_Op(cinq,nloop);
  Bench_10_Op(six,nloop);
  Bench_10_Op(sept,nloop);
  Bench_10_Op(huit,nloop);
  Bench_10_Op(neuf,nloop);
  Bench_10_Op(dix,nloop);
  Bench_10_Op(onze,nloop);
  Bench_10_Op(douze,nloop);
  Bench_10_Op(treize,nloop);
  Bench_10_Op(quatorze,nloop);
  Bench_10_Op(quinze,nloop);
  Bench_10_Op(seize,nloop);
  Bench_10_Op(dixsept,nloop);
  Bench_10_Op(dixhuit,nloop);
  Bench_10_Op(dixneuf,nloop);
  Bench_10_Op(vingt,nloop);

  Bench_20_Op(un,nloop);
  Bench_20_Op(deux,nloop);
  Bench_20_Op(trois,nloop);
  Bench_20_Op(quatre,nloop);
  Bench_20_Op(cinq,nloop);
  Bench_20_Op(six,nloop);
  Bench_20_Op(sept,nloop);
  Bench_20_Op(huit,nloop);
  Bench_20_Op(neuf,nloop);
  Bench_20_Op(dix,nloop);
  Bench_20_Op(onze,nloop);
  Bench_20_Op(douze,nloop);
  Bench_20_Op(treize,nloop);
  Bench_20_Op(quatorze,nloop);
  Bench_20_Op(quinze,nloop);
  Bench_20_Op(seize,nloop);
  Bench_20_Op(dixsept,nloop);
  Bench_20_Op(dixhuit,nloop);
  Bench_20_Op(dixneuf,nloop);
  Bench_20_Op(vingt,nloop);


  return 0;
}


template <int Num> void Bench_1_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  for (j=0; j<ntest; j++){

    TFad< Num > x1(2.f,0), x2(3.f), x4;

     Timer mytime;
     mytime.start();

     for (i=0; i<nloop; ++i) 
       x4 = x1 + x2;

     mytime.stop();
     time1 += mytime.elapsedSeconds();
    //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


    TinyFad< Num > y1(2.,0), y2((float)3.), y4;

    mytime.start();

    for (i=0; i<nloop; ++i)
       y4 = y1 + y2;

    mytime.stop();
    time2 += mytime.elapsedSeconds();
    //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;


    float z1=2., z2=3., z4;
    float tz1[Num], tz2[Num], tz4[Num];

    mytime.start();

    for (i=0; i<nloop; ++i) {
      z4 = z1 + z2;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    mytime.start();

    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val() + y2.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k) + y2.dx(k);
    }

    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;


  }

  cout.width(12);
  cout << (time4 / time3);
  cout.width(12);
  cout << (time1 / time3) ;
  cout.width(12);
  cout << (time2 / time3) << "\n";

}


template <int Num> void Bench_2_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


   TFad< Num > x1(2.,0), x2(3.), x3(5.), x4, x5(3.);

  for (j=0; j<ntest; j++){

     Timer mytime;
     mytime.start();

     for (i=0; i<nloop; ++i) 
       x4 = x1+x2+x3;

     mytime.stop();
     time1 += mytime.elapsedSeconds();
     //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


    TinyFad< Num > y1(2.,0), y2((float)3.), y3((float)5.), y4, y5((float)3.);

    mytime.start();

    for (i=0; i<nloop; ++i)
       y4 = y1+y2+y3;

    mytime.stop();
    time2 += mytime.elapsedSeconds();
    //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;


    float z1=2., z2=3., z3=5., z4, z5=3.;
    float tz1[Num], tz2[Num], tz3[Num], tz4[Num], tz5[Num];

    mytime.start();

    for (i=0; i<nloop; ++i) {
      z4 = z1+z2+z3;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k]+tz3[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    mytime.start();

    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val()+y2.val()+y3.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k)+y2.dx(k)+y3.dx(k);
    }

    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;


  }

  cout.width(12);
  cout << (time4 / time3);
  cout.width(12);
  cout << (time1 / time3) ;
  cout.width(12);
  cout << (time2 / time3) << "\n";

}


template <int Num> void Bench_3_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  TFad<Num > x1(2.,0), x2(3.), x3(5.), x4, x5(3.);

  for (j=0; j<ntest; j++){

     Timer mytime;
     mytime.start();

     for (i=0; i<nloop; ++i) 
       x4 = x1+x2+x3+x3;

     mytime.stop();
     time1 += mytime.elapsedSeconds();
     //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


    TinyFad< Num > y1(2.,0), y2((float)3.), y3((float)5.), y4, y5((float)3.);

    mytime.start();

    for (i=0; i<nloop; ++i)
       y4 = y1+y2+y3+y3;

    mytime.stop();
    time2 += mytime.elapsedSeconds();
    //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;


    float z1=2., z2=3., z3=5., z4, z5=3.;
    float tz1[Num], tz2[Num], tz3[Num], tz4[Num], tz5[Num];

    mytime.start();

    for (i=0; i<nloop; ++i) {
      z4 = z1+z2+z3+z3;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k]+tz3[k]+tz3[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    mytime.start();

    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val()+y2.val()+y3.val()+y3.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k);
    }

    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;


  }

  cout.width(12);
  cout << (time4 / time3);
  cout.width(12);
  cout << (time1 / time3) ;
  cout.width(12);
  cout << (time2 / time3) << "\n";

}


template <int Num> void Bench_4_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  TFad< Num > x1(2.,0), x2(3.), x3(5.), x4, x5(3.);

  for (j=0; j<ntest; j++){

    Timer mytime;
    mytime.start();
    
    for (i=0; i<nloop; ++i) 
      x4 = x1+x2+x3+x3+x5;
    
    mytime.stop();
    time1 += mytime.elapsedSeconds();
    //     //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


    TinyFad< Num > y1(2.f,0), y2(3.f), y3(5.f), y4, y5(3.f);
    mytime.start();
    
    for (i=0; i<nloop; ++i)
      y4 = y1+y2+y3+y3+y5;
    
    mytime.stop();
    time2 += mytime.elapsedSeconds();
    //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;


    float z1=2., z2=3., z3=5., z4, z5=3.;
    float tz1[Num], tz2[Num], tz3[Num], tz4[Num], tz5[Num];
    mytime.start();
    
    for (i=0; i<nloop; ++i) {
      z4 = z1+z2+z3+z3+z5;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    
    mytime.start();
    
    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val()+y2.val()+y3.val()+y3.val()+y5.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k);
    }
    
    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;
    
    
  }

  cout << setw(12) << time1 << setw(12) << time2 << setw(12) << time3 << setw(12) << time4 << endl;

}


template <int Num> void Bench_10_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  TFad< Num > x1(2.,0), x2(3.), x3(5.), x4, x5(3.);

  for (j=0; j<ntest; j++){

     Timer mytime;
     mytime.start();

     for (i=0; i<nloop; ++i) 
       x4 = x1+x2+x3+x3+x5+x1+x2+x3+x3+x5;

     mytime.stop();
     time1 += mytime.elapsedSeconds();
     //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


     TinyFad< Num > y1(2.,0), y2((float)3.), y3((float)5.), y4, y5((float)3.);

     mytime.start();

     for (i=0; i<nloop; ++i)
       y4 = y1+y2+y3+y3+y5+y1+y2+y3+y3+y5;

     mytime.stop();
     time2 += mytime.elapsedSeconds();
     //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;

     
    float z1=2., z2=3., z3=5., z4, z5=3.;
    float tz1[Num], tz2[Num], tz3[Num], tz4[Num], tz5[Num];

    mytime.start();

    for (i=0; i<nloop; ++i) {
      z4 = z1+z2+z3+z3+z5+z1+z2+z3+z3+z5;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k]+tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    mytime.start();

    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val()+y2.val()+y3.val()+y3.val()+y5.val()+y1.val()+y2.val()+y3.val()+y3.val()+y5.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k)+y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k);
    }

    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;


  }

  cout << setw(12) << time1 << setw(12) << time2 << setw(12) << time3 << setw(12) << time4 << endl;

}

template <int Num> void Bench_20_Op(const TinyFad<Num> & x, const int nloop){
  int i,j,k;
  const int ntest = 3;

  cout.setf(ios::fixed,ios::floatfield);
  cout.width(4);
  cout << Num;
  cout.flush();

  float time1 = 0., time2 = 0., time3 = 0., time4 = 0.;


  TFad< Num > x1(2.,0), x2(3.), x3(5.), x4, x5(3.);

  for (j=0; j<ntest; j++){

    Timer mytime;
    mytime.start();
     
    for (i=0; i<nloop; ++i) 
      x4 = x1+x2+x3+x3+x5+x1+x2+x3+x3+x5+x1+x2+x3+x3+x5+x1+x2+x3+x3+x5;

    mytime.stop();
    time1 += mytime.elapsedSeconds();
    //cout << "TFad : " << mytime.elapsedSeconds();// << endl;


    TinyFad< Num > y1(2.,0), y2((float)3.), y3((float)5.), y4, y5((float)3.);

    mytime.start();

    for (i=0; i<nloop; ++i)
       y4 = y1+y2+y3+y3+y5+y1+y2+y3+y3+y5+y1+y2+y3+y3+y5+y1+y2+y3+y3+y5;

    mytime.stop();
    time2 += mytime.elapsedSeconds();
    //cout << " , TinyFad : " << mytime.elapsedSeconds();// << endl;


    float z1=2., z2=3., z3=5., z4, z5=3.;
    float tz1[Num], tz2[Num], tz3[Num], tz4[Num], tz5[Num];

    mytime.start();

    for (i=0; i<nloop; ++i) {
      z4 = z1+z2+z3+z3+z5+z1+z2+z3+z3+z5+z1+z2+z3+z3+z5+z1+z2+z3+z3+z5;
      for (k=0; k<Num; ++k)
	tz4[k] = tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k]+tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k]+tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k]+tz1[k]+tz2[k]+tz3[k]+tz3[k]+tz5[k];
    }

    mytime.stop();
    time3 += mytime.elapsedSeconds();
    //cout << " , float : " << mytime.elapsedSeconds();// << endl;

    mytime.start();

    for (i=0; i<nloop; ++i) {
      y4.val() = y1.val()+y2.val()+y3.val()+y3.val()+y5.val()+y1.val()+y2.val()+y3.val()+y3.val()+y5.val()+y1.val()+y2.val()+y3.val()+y3.val()+y5.val()+y1.val()+y2.val()+y3.val()+y3.val()+y5.val();
      for (k=0; k<Num; ++k)
	y4.dx(k) = y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k)+y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k)+y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k)+y1.dx(k)+y2.dx(k)+y3.dx(k)+y3.dx(k)+y5.dx(k);
    }

    mytime.stop();
    time4 += mytime.elapsedSeconds();
    //cout << " , Acess Members : " << mytime.elapsedSeconds() << endl;


  }

  cout << setw(12) << time1 << setw(12) << time2 << setw(12) << time3 << setw(12) << time4 << endl;

}

