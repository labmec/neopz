#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

double pi(double *u,double *v);
void pv(double *u,double *v,double *w);

/*        !! ATENÇÃO !!
   Este programa faz duas operações:
   1.- produto interno entre u e v (pi)
   2.- produto vetorial entre u e v (pv)
   continua=0 procesa dois vetores
   continua=1 lê mais um vetor
   continua=2 lê mais dois pontos
   o programa lê: u,v ou u,v,w ou u,v,w,p
   continua=0 faz o produto entre u e v
   continua=1 faz o produto entre v-u e w-u
   continua=2 faz o produto entre v-u e p-w
*/

void main(){

    ifstream vec("vectors.txt");
    double u[3],v[3],w[3],p[3];
    int continua,instante;
    do{
       w[0]=0.; w[1]=0.; w[2]=0.;
       vec >> u[0] >> u[1] >> u[2];
       vec >> v[0] >> v[1] >> v[2];
       vec >> continua;
       if(continua==1){
          vec >> w[0] >> w[1] >> w[2];
       }else if(continua==2){
          vec >> w[0] >> w[1] >> w[2];
          vec >> p[0] >> p[1] >> p[2];
       }
       if(continua==1){
          double a[3];
          a[0] = u[0]; a[1] = u[1]; a[2] = u[2];
          u[0] = v[0] - a[0];
          u[1] = v[1] - a[1];
          u[2] = v[2] - a[2];
          v[0] = w[0] - a[0];
          v[1] = w[1] - a[1];
          v[2] = w[2] - a[2];
          cout << "\nVetores : \n";
          cout << "u = ( "  << u[0] << " , " << u[1]
               << " , " << u[2] << " )" << endl;
          cout << "v = ( "  << v[0] << " , " << v[1]
               << " , " << v[2] << " )" << endl;
       } else if(continua==2){
          u[0] = u[0] - v[0];
          u[1] = u[1] - v[1];
          u[2] = u[2] - v[2];
          v[0] = w[0] - p[0];
          v[1] = w[1] - p[1];
          v[2] = w[2] - p[2];
          cout << "\nVetores : \n";
          cout << "u = ( "  << u[0] << " , " << u[1]
               << " , " << u[2] << " )" << endl;
          cout << "v = ( "  << v[0] << " , " << v[1]
               << " , " << v[2] << " )" << endl;
       }
       double r = pi(u,v);
       pv(u,v,w);
       //ofstream result("result.txt");
       //result << r;
       cout << "\nproduto interno : " << r;
       cout << "\nproduto vetorial : " << "( " << w[0]
            << " , " << w[1] << " , " << w[2] << " )" << endl;
       cin >> instante;
    } while(instante);
}

double pi(double *u,double *v){
     return (u[0]*v[0]+u[1]*v[1]+u[2]*v[2]);
}

void pv(double *u,double *v,double *w){

    w[0] = u[1]*v[2] - u[2]*v[1];
    w[1] = u[2]*v[0] - u[0]*v[2];
    w[2] = u[0]*v[1] - u[1]*v[0];
}

