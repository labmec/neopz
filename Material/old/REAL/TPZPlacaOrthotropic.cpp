/**
 * \file
 * @brief Contains implementations of the TPZPlacaOrthotropic methods.
 */
#include "pzintel.h" 
#include "TPZPlacaOrthotropic.h"
#include "pzgnode.h"
#include "pzgeoel.h"
#include "pzquad.h"
//#include "pztempmat.h"
#include "TPZMaterial.h"

using namespace std;

TPZPlacaOrthotropic::TPZPlacaOrthotropic(TPZGeoEl *gel,REAL zmin, REAL zmax ){ 
	
	fGeoEl = gel;
	fH = zmax-zmin;
	fZMin = zmin;
	fZMax = zmax;
	fTensorVar = -1;
	fIntel = dynamic_cast<TPZInterpolatedElement *>(gel->Reference());
	if(fIntel) {
		fTensorVar = fIntel->Material()->VariableIndex("Tensor");
	}
}

TPZPlacaOrthotropic::TPZPlacaOrthotropic(){ 
	
	fGeoEl = 0;
	fH = 0.;
	fZMin = 0.;
	fZMax = 0.;
	fIntel = 0;
	fTensorVar = -1;
}

void TPZPlacaOrthotropic::Tensor(TPZVec<REAL> &ksi, TPZFMatrix<REAL> &T) {
	if(fTensorVar == -1) {
		if(!fIntel) fIntel = dynamic_cast<TPZInterpolatedElement *>(fGeoEl->Reference());
		if(!fIntel) return;
		fTensorVar = fIntel->Material()->VariableIndex("Tensor");
	}
	if(fTensorVar == -1) return;
	//  REAL ksi = -1+2.*(z-fZMin)/(fZMax-fZMin);
	//  if(ksi < -1. || ksi > 1.) return;
	//  TPZManVector<REAL,3> co(3,0.);
	//  co[2] = ksi;
	TPZManVector<STATE> tensor(9);
	fIntel->Solution(ksi,fTensorVar,tensor);
	int i;
	for(i=0; i<9; i++) {//original
		T(i%3,i/3) = tensor[i];
	}
	//   for(i=0; i<3; i++) {//�sim�rico
	//     for(j=0; j<3; j++) {
	//       T(i,j) = tensor[3*i+j];
	//     }
	//   }
}

REAL TPZPlacaOrthotropic::Moment(REAL zref, TPZVec<REAL> &normal, TPZVec<REAL> &direction){
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZManVector<REAL,3> pos(1,0.), ksi(3,0.);
	TPZFNMatrix<9> tensor(3,3,0.);
	int ip;
	REAL moment = 0.;
	for(ip = 0; ip<np; ip++) {
		REAL weight;
		rule.Point(ip,pos,weight);
		ksi[2] = pos[0];
		Tensor(ksi,tensor);
		REAL z = fZMin + (fZMax-fZMin)*(pos[0]+1.)/2.;//original
		//REAL f = (fZMin + fZMax)/2.0;//n�
		//REAL z = fH*pos[0]/2. + f;//n�
		REAL tension = 0.;
		int n,d;
		for(n=0; n<3; n++) {
			for(d=0; d<3; d++) {
				tension += normal[n]*tensor(n,d)*direction[d];
			}
		}
		moment += weight*tension*(z-zref);//original
		//moment += weight*tension*z;//n�
	}
	moment *= fH/2.;
	return moment;
	
}

REAL TPZPlacaOrthotropic::Force(TPZVec<REAL> &normal, TPZVec<REAL> &direction){
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZManVector<REAL,3> pos(1,0.), ksi(3,0.);
	TPZFNMatrix<9> tensor(3,3,0.);
	int ip;
	REAL force = 0.;
	for(ip = 0; ip<np; ip++) {
		REAL weight;
		rule.Point(ip,pos,weight);
		ksi[2] = pos[0];
		Tensor(ksi,tensor);
		REAL tension = 0.;
		int n,d;
		for(n=0; n<3; n++) {
			for(d=0; d<3; d++) {
				tension += normal[n]*tensor(n,d)*direction[d];
			}
		}
		force += weight*tension;
	}
	force *= fH/2.;
	return force;
}

void TPZPlacaOrthotropic::Print(){
	
	fIntel->Print();
	cout << "Plate height : " ;
	cout << this->Height();//this �o objeto placa
	
}

//se fosse TPZMulticamadaOrthotropic::Print, o this seria multcam
//se a func� fosse declarada como void Print(){...}, n� haveria This.
//caso chamemos a func� Print com o objeto placa2, ao inv� de placa->Print, o This �o objeto placa2.


void TPZPlacaOrthotropic::IdentifyCompEl() {
	fIntel = dynamic_cast<TPZInterpolatedElement *>(fGeoEl->Reference());
	if(fIntel) {
		fTensorVar = fIntel->Material()->VariableIndex("Tensor");
	}
	
}

REAL TensionNorm(TPZFMatrix<REAL> &tension,int dimrow,int dimcol);
void TPZPlacaOrthotropic::PrintTensors(std::ostream &out,TPZFMatrix<REAL> &tensorin,TPZFMatrix<REAL> &tensorout) {
	
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZFMatrix<REAL> TensorOut(np,9),DiffTension(np,9);
	TPZManVector<REAL,3> pos(1,0.), ksi(3,0.), x(3,0.);
	TPZFNMatrix<9> tensor(3,3,0.);
	int ip,i,j;
	for(ip = 0; ip<np; ip++) {
		REAL weight;
		rule.Point(ip,pos,weight);
		ksi[2] = pos[0];
		fGeoEl->X(ksi,x);
		Tensor(ksi,tensor);
		out << "Tensor at pos " << pos[0] << " x " << x << endl;
		for(i=0; i<3; i++) for(j=0; j<3; j++) tensorout(ip,3*i+j) = tensor(i,j);
		for(j=0; j<9; j++) DiffTension(ip,j) = tensorout(ip,j) - tensorin(ip,j);
	}
	out << " tensor(n) - tensor(n-1) = " << endl;
	for(i=0; i<np; i++){
		for(j=0; j<9; j++){
			out <<  DiffTension(i,j) << " ";
		}
		out << endl;
	}
	REAL normat = TensionNorm(DiffTension,np,9);
	out << "Norma of tensor out - tensor in = " << normat << endl;
}

void TPZPlacaOrthotropic::PrintTensors(std::ostream &out) {
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZManVector<REAL,3> pos(1,0.), ksi(3,0.), x(3,0.);
	TPZFNMatrix<9> tensor(3,3,0.);
	int ip;
	REAL weight;
	TPZStack<REAL> point;
	point.Push(-1.0);
	for(ip = 0; ip<np; ip++){
		rule.Point(ip,pos,weight);
		point.Push(pos[0]);
	}
	point.Push(1.0);
	np += 2;
	for(ip = 0; ip<np; ip++) {
		
		ksi[2] = point[ip];
		fGeoEl->X(ksi,x);
		Tensor(ksi,tensor);
		out << "Tensor at pos " << point[ip] << " x " << x << " -> ";
		int i,j;
		for(i=0; i<3; i++) {
			for(j=0; j<3; j++){
				REAL tension = tensor(i,j);
				if(fabs(tension) < 1.e-10) out << 0 << " ";
				else out << tension << " ";
			}
			//out << endl;
		}
		out << endl;
		REAL normat = TensionNorm(tensor,3,3);
		//    out << "Norma do tensor " << normat << endl;
		cout << "Norma do tensor " << normat << endl;
	}
}

/** @brief Returns norm of the tension */
REAL TensionNorm(TPZFMatrix<REAL> &tension,int dimrow,int dimcol) {
	
	int i,j;
	REAL val = 0.0;
	for(i=0;i<dimrow;i++) for(j=0;j<dimcol;j++) val += tension(i,j)*tension(i,j);
	return sqrt(val);
}

REAL TPZPlacaOrthotropic::GradMoment(REAL zref, TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction){
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZManVector<REAL,3> pos(1,0.), ksi1(3,0.), ksi2(3,0.), x1(3,0.), x2(3,0.);
	TPZFNMatrix<9> tensor1(3,3,0.), tensor2(3,3,0.), tensordif(3,3,0.);
	ksi1[0] = -graddir[0];
	ksi1[1] = -graddir[1];
	ksi2[0] = graddir[0];
	ksi2[1] = graddir[1];
	fGeoEl->X(ksi1,x1);
	fGeoEl->X(ksi2,x2);
	REAL dist = sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]));
	int ip;
	REAL moment = 0.;
	for(ip = 0; ip<np; ip++) {
		REAL weight;
		rule.Point(ip,pos,weight);
		ksi1[2] = pos[0];
		ksi2[2] = pos[0];
		Tensor(ksi1,tensor1);
		Tensor(ksi2,tensor2);
		tensordif = tensor2;
		tensordif -= tensor1;
		REAL z = fZMin + (fZMax-fZMin)*(pos[0]+1.)/2.;// = (z+f)
		REAL tension = 0.;
		int n,d;
		for(n=0; n<3; n++) {
			for(d=0; d<3; d++) {
				tension += normal[n]*tensordif(n,d)*direction[d];
			}
		}
		moment += weight*tension*(z-zref);//original
		//moment += weight*tension*z;//n�
	}
	moment *= fH/(2.*dist);
	return moment;
	
}

REAL TPZPlacaOrthotropic::GradForce(TPZVec<REAL> &graddir, TPZVec<REAL> &normal, TPZVec<REAL> &direction){
	TPZInt1d rule(8);
	int np = rule.NPoints();
	TPZManVector<REAL,3> pos(1,0.), ksi1(3,0.), ksi2(3,0.), x1(3,0.), x2(3,0.);
	TPZFNMatrix<9> tensor1(3,3,0.), tensor2(3,3,0.), tensordif(3,3,0.);
	ksi1[0] = -graddir[0];
	ksi1[1] = -graddir[1];
	ksi2[0] = graddir[0];
	ksi2[1] = graddir[1];
	fGeoEl->X(ksi1,x1);
	fGeoEl->X(ksi2,x2);
	REAL dist = sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]));
	int ip;
	REAL force = 0.;
	for(ip = 0; ip<np; ip++) {
		REAL weight;
		rule.Point(ip,pos,weight);
		ksi1[2] = pos[0];
		ksi2[2] = pos[0];
		Tensor(ksi1,tensor1);
		Tensor(ksi2,tensor2);
		tensordif = tensor2;
		tensordif -= tensor1;
		REAL tension = 0.;
		int n,d;
		for(n=0; n<3; n++) {
			for(d=0; d<3; d++) {
				tension += normal[n]*tensordif(n,d)*direction[d];
			}
		}
		force += weight*tension;
	}
	force *= fH/(2.*dist);
	return force;
}

void TPZPlacaOrthotropic::GradTensor(TPZVec<REAL> &graddir, TPZVec<REAL> &ksi, TPZFMatrix<REAL> &T) {
	if(fTensorVar == -1) {
		if(!fIntel) fIntel = dynamic_cast<TPZInterpolatedElement *>(fGeoEl->Reference());
		if(!fIntel) return;
		fTensorVar = fIntel->Material()->VariableIndex("Tensor");
	}
	if(fTensorVar == -1) return;
	//  REAL ksi = -1+2.*(z-fZMin)/(fZMax-fZMin);
	//  if(ksi < -1. || ksi > 1.) return;
	//  TPZManVector<REAL,3> co(3,0.);
	//  co[2] = ksi;
	TPZManVector<REAL,3> pos(1,0.), ksi1(3,0.), ksi2(3,0.), x1(3,0.), x2(3,0.);
	TPZManVector<STATE> tensor1(0,0.), tensor2(9,0.);
	ksi1[0] = ksi[0]-graddir[0];
	ksi1[1] = ksi[1]-graddir[1];
	ksi1[2] = ksi[2]-graddir[2];
	ksi2[0] = ksi[0]+graddir[0];
	ksi2[1] = ksi[1]+graddir[1];
	ksi2[2] = ksi[2]+graddir[2];
	fGeoEl->X(ksi1,x1);
	fGeoEl->X(ksi2,x2);
	REAL dist = sqrt((x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]));
	TPZManVector<REAL> tensor(9);
	fIntel->Solution(ksi1,fTensorVar,tensor1);
	fIntel->Solution(ksi2,fTensorVar,tensor2);
	int i;
	for(i=0; i<9; i++) {//�sim�rico OK
		T(i%3,i/3) = (tensor2[i]-tensor1[i])/dist;
	}
	//   for(i=0; i<3; i++) {//n� sim�rico
	//     for(j=0; j<3; j++) {
	//       T(i,j) = (tensor2[3*i+j]-tensor1[3*i+j])/dist;
	//     }
	//   }
}
