/// ---- GLOBALS ----

#define NTHREADS 2
#define MAXLVL 5
#define NX 5
#define NY 5

/// -----------

//#define MR

#define REFTHREADV3 				//main refthread third version
//#define REFTHREADV2				//main refthread second version
//#define REFTHREADV1				//main refthread first working version
//#define REFTHREADSERIAL           //serial implementation
//#define ADMCHUNKTESTS				//enabling the tests on admchunkvectorthreadsafe

/// -----------

/// On ADMCHUNKTESTS:

	#define TEST 7

	// if TEST = 1				//for testing operator=() on the vector
	// if TEST = 2				//for testing AllocateNewElement() on the vector
	// if TEST = 3				//for testing Resize() on the vector
	// if TEST = 4				//for testing SetFree() on the vector
	// if TEST = 5				//for testing operator[]() on the vector
	// if TEST = 6				//for testing FindObject() on the vector
	// if TEST = 7				//for testing CompactDataStructure() on the vector

/// -----------

/// ---- GLOBALS ----

#include <cmath>
#include <iostream>
#include <pz_gettime.h>
#include <sstream>
#include <thread>
#include <pzgmesh.h>
#include <pzgeoel.h>
#include "pzvec.h"
#include "TPZVTKGeoMesh.h"
#include "TPZRefPatternDataBase.h"

#include "pzskylstrmatrix.h"
#include "pzsolve.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "pzgeoelbc.h"

using namespace std;

double get_time()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	double d = t.tv_sec + (double) t.tv_usec/1000000;
	return d;
}	// Simple function for taking time

TPZGeoMesh* GetMesh (int nx, int ny)
{
	int i,j;
	int64_t id, index;
	
	REAL lx = 1.;
	REAL ly = 1.;
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	TPZVec <REAL> coord (3,0.);
	
	for(i = 0; i < nx; i++){
		for(j = 0; j < ny; j++){
			id = i*ny + j;
			coord[0] = (i)*lx/(nx - 1);
			coord[1] = (j)*ly/(ny - 1);
			coord[2] = 0.;
			index = gmesh->NodeVec().AllocateNewElement();
			gmesh->NodeVec()[index] = TPZGeoNode(id,coord,*gmesh);
		}
	}
	
	TPZGeoEl * elvec[(const int)((nx-1)*(ny-1))];
	
	TPZVec <int64_t> connect(4,0);
	
	for(i = 0; i < (nx - 1); i++){
		for(j = 0; j < (ny - 1); j++){
			index = (i)*(ny - 1)+ (j);
			connect[0] = (i)*ny + (j);
			connect[1] = connect[0]+(ny);
			connect[2] = connect[1]+1;
			connect[3] = connect[0]+1;
			elvec[index] = gmesh->CreateGeoElement(EQuadrilateral,connect,1,id);
		}
	}
	
	gmesh->BuildConnectivity();
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	
	return gmesh;
}

#ifdef MR

TPZVec <TPZGeoEl *> sons;
double err;

int option, maxlvl, nx, ny, nthreads;

struct divide_data_parallel {
	int el;
	TPZGeoMesh* mesh;
};

void *divide_parallel (void *arg)
{
	divide_data_parallel *idata = (divide_data_parallel*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	delete idata;
	
	TPZVec <TPZGeoEl *> sons;
    
    int iref, nref = maxlvl;
	for (iref=0; iref<nref; iref++)
	{
		int nelements = imesh->NElements();
		
		for (int iel=0; iel<nelements; iel++)
		{
			TPZGeoEl * gel = imesh->ElementVec()[iel];
			
#ifdef PZDEBUG
			if(!gel)
			{
				DebugStop();
			}
#endif
			
			if(!gel->HasSubElement())
			{
				gel->Divide(sons);
			}
		}
	}
    
}

int main_parallel ()
{	
    std::vector<std::thread> threads;
    
    TPZGeoMesh *gmeshes[nthreads];
    for (int i=0; i<nthreads; i++) {
        gmeshes[i] = GetMesh(nx, ny);
    }
	
    int err;
    
	double time_start = get_time();
	//cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS (MRPARALLEL)***\n";
	
    for (int i=0; i<nthreads; i++) {
        divide_data_parallel *sdata;
        sdata = new divide_data_parallel;
        sdata->mesh = gmeshes[i];
        threads.push_back(std::thread(divide_parallel,sdata));
    }
    
    for (int i=0; i<nthreads; i++) {
      threads[i].join();        
    }
    
	double time_end = get_time();
    //cout << "\n****END****\n";
	//cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
    
    for (int i=0; i<nthreads; i++) { 
        string filename;
        
        ostringstream ss;
        ss << i;
        string filenumber = ss.str();
       
        ostringstream ss2;
        ss2 << nthreads;
        string nthreads_str = ss2.str();
        
        filename = "after_PTHREAD_";
        filename.insert(filename.size(), nthreads_str);
        filename.insert(filename.size(),"meshes_mesh_");
        filename.insert(filename.size(), filenumber);
        filename.insert(filename.size(),".vtk");
        
        ofstream af(filename.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshes[i], af);
        af.close();
    }
    
    ofstream af;
    af.open("after_PTHREAD_report.txt", ios::app);
    
    af << option << " " << maxlvl << " " << nx << " " << ny << " " << nthreads << " " << (time_end-time_start) << "\n";
    
    af.close();
    
	return 0;
}

int main_serial ()
{ 	
	TPZGeoMesh *gmeshes[nthreads];
    for (int i=0; i<nthreads; i++) {
        gmeshes[i] = GetMesh(nx, ny);
    }
	
	double time_start = get_time();
    //cout << "\n\n***STARTING SERIAL REFINEMENT PROCESS (MRSERIAL)***\n";
    
    for (int i=0; i<nthreads; i++) {
        int iref, nref = maxlvl;
        for (iref=0; iref<nref; iref++)
        {
            int nelements = gmeshes[i]->NElements();
            
            for (int iel=0; iel<nelements; iel++)
            {
                TPZGeoEl * gel = gmeshes[i]->ElementVec()[iel];
                
#ifdef PZDEBUG
                if(!gel)
                {
                    DebugStop();
                }
#endif
                
                if(!gel->HasSubElement())
                {
                    gel->Divide(sons);
                }
            }
        }
	}
	
	double time_end = get_time();
    //cout << "\n****END****\n";
	//cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
	    
    for (int i=0; i<nthreads; i++) { 
        string filename;
        
        ostringstream ss;
        ss << i;
        string filenumber = ss.str();
        
        ostringstream ss2;
        ss2 << nthreads;
        string nthreads_str = ss2.str();
        
        filename = "after_SERIAL_";
        filename.insert(filename.size(), nthreads_str);
        filename.insert(filename.size(),"meshes_mesh_");
        filename.insert(filename.size(), filenumber);
        filename.insert(filename.size(),".vtk");
        
        ofstream af(filename.c_str());
        TPZVTKGeoMesh::PrintGMeshVTK(gmeshes[i], af);
        af.close();
    }

    ofstream af;
    af.open("after_SERIAL_report.txt", ios::app);
    
    af << option << " " << maxlvl << " " << nx << " " << ny << " " << nthreads << " " << (time_end-time_start) << "\n";
    
    af.close();
    
	return 0;
}

int main (int argc, char *argv[]) {
    
    if (argc!=6) {
        cout << "\n ******* Error! Invalid parameters! *******\n";
        return 0;
    }
    
    option = atoi(argv[1]);
    
    maxlvl = atoi(argv[2]);
    
    nx = atoi(argv[3]);
    
    ny = atoi(argv[4]);
    
    nthreads = atoi(argv[5]);
    
    if (option==1) {
        main_serial();
    }
    
    else if (option==2) {
        main_parallel();
    }
}

#endif

//------------------------------------------------------------//

#ifdef REFTHREADV3

#include <queue>
#include <mutex>
std::vector<std::thread> threads;
std::mutex idslock;

TPZVec <int> *ids;
queue <int> *todo[NTHREADS];

double err;
const int azero = 0;

struct divide_data {
	int threadnumber;
	TPZGeoMesh* mesh;
};

TPZStack <int> getneighbours (TPZGeoMesh* mesh, int iel) 
{
	//double getneigh_time_start = get_time();
	TPZStack <int> nbs;
	int countneighbours=0, nneighbours=0;
	
	int nsides = mesh->ElementVec()[iel]->NSides();
	for (int i=0; i<nsides; i++)
	{
		TPZGeoElSide neighbour = mesh->ElementVec()[iel]->Neighbour(i);
		TPZGeoElSide thisside (mesh->ElementVec()[iel],i);
		if (neighbour.Exists()) 
		{
			int count=0;
			while (neighbour!=thisside&&count++<30)
			{
				countneighbours = 0;
				nneighbours = nbs.size();
				for (int j=0; j<nneighbours; j++)
				{
					if (neighbour.Element()->Id()!=nbs[j]) countneighbours++;
				}
				if (countneighbours==nneighbours)
				{
					int pos = neighbour.Element()->Id();
					nbs.Push(pos);
				}
				neighbour = neighbour.Neighbour();
			}
		}
	}
	//double getneigh_time_end = get_time();
	//cout << "getneighbours() took "<< (getneigh_time_end-getneigh_time_start) << " seconds.\n\n";
    
	return nbs;
}

int lockneighbours (TPZGeoMesh *mesh, int iel)
{
	//double lockneigh_time_start = get_time();
	std::scoped_lock<std::mutex> lck(idslock);
	TPZStack <int> neighbours = getneighbours(mesh, iel);
	
	while (neighbours.size()) {
		int pos = neighbours.Pop();
		ids->Fill((int) 1, pos, 1);
	}
	
	ids->Fill((int) 1, iel, 1);
	
	//double lockneigh_time_end = get_time();
	//cout << "lockneighbours() took "<< (lockneigh_time_end-lockneigh_time_start) << " seconds.\n\n";
	
	return 0;
}

int unlockneighbours (TPZGeoMesh *mesh, int iel, int threadnumber, TPZVec <TPZGeoEl*> sons)
{
  std::scoped_lock lck(idslock);
	
	TPZStack <int> neighbours = getneighbours(mesh, iel);
	
	ids->Fill(0, iel, 1);
	
	while (neighbours.size()) {
		int pos = neighbours.Pop();
		ids->Fill(0, pos, 1);
	} 
	
	
    for (int i=0; i<sons.NElements(); i++) {
		if (sons[i]->Level()<MAXLVL) {
			todo[threadnumber]->push(sons[i]->Id());
		}
        else break;
	}
    
	//double unlockneigh_time_end = get_time();
	//cout << "unlockneighbours() took "<< (unlockneigh_time_end-unlockneigh_time_start) << " seconds.\n\n";
	
	return 0;
}

bool checkneighbours (TPZGeoMesh* mesh, int iel) 
{
	//double checkneigh_time_start = get_time();
	
    if (ids->operator[](iel))
	{
		//double checkneigh_time_end = get_time();
		//cout << "checkneighbours("<< iel<< ")=true took "<< (checkneigh_time_end-checkneigh_time_start) << " seconds.\n\n";
		return true;
	}
	
	else
	{
		int countneighbours=0;
		int nneighbours=0;
		
		TPZStack <int> neighbours = getneighbours(mesh, iel);
		
		nneighbours = neighbours.size();
		
		while (neighbours.size())
		{
			int neigh = neighbours.Pop();
			if(!ids->operator[](neigh)) countneighbours++;
		}
		
		if (countneighbours==nneighbours)
        {
			//double checkneigh_time_end = get_time();
			//cout << "checkneighbours("<< iel<< ")=false took "<< (checkneigh_time_end-checkneigh_time_start) << " seconds.\n\n";
			return false;
		}
		
		else
		{
			//double checkneigh_time_end = get_time();
			//cout << "checkneighbours("<< iel<< ")=true took "<< (checkneigh_time_end-checkneigh_time_start) << " seconds.\n\n";
			return true;
		}
	}
}

void *divide(void *arg)
{
    divide_data *idata = (divide_data*) arg;
    TPZGeoMesh *imesh = idata->mesh;
    int threadnumber = idata->threadnumber;
    delete idata;
    std::unique_lock lck(idslock);
    TPZVec <TPZGeoEl *> sons;
    
    while (todo[threadnumber]->size())
      {
        int iel = todo[threadnumber]->front();
        todo[threadnumber]->pop();
        
        lck.lock();
        
        if (!checkneighbours(imesh, iel))
          {
            lockneighbours(imesh, iel);
            
            imesh->ElementVec()[iel]->Divide(sons);
            
            unlockneighbours(imesh, iel, threadnumber, sons);
          }
        
        else {
          todo[threadnumber]->push(iel);
          lck.unlock();
        }
      }
    lck.unlock();
}

int main ()
{	
	TPZGeoMesh *gmesh = GetMesh(NX, NY);
	
	//ofstream bef("before_PTHREAD3.vtk");
	//TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file
    
	int threadnumber = 0;
  int nelementsperthread;
	REAL totalElements = 1;
	REAL alevel;
	
	// start of initialization
	for (alevel=1; alevel<=MAXLVL; alevel++)
	{
		totalElements+=pow((REAL)4., alevel);
	}
	
	totalElements*=gmesh->NElements();
	
	ids = new TPZVec <int> (totalElements, azero);
    
    for (int i=0; i<NTHREADS; i++)
    {
        todo[i] = new queue <int>;
    }
    
    nelementsperthread = gmesh->NElements()/NTHREADS;
	
    for (threadnumber=0; threadnumber<NTHREADS; threadnumber++)
    {
        for (int i=(nelementsperthread*threadnumber); i<(nelementsperthread*(threadnumber+1)); i++)
        {
            todo[threadnumber]->push(gmesh->ElementVec()[i]->Id());
            if (threadnumber==(NTHREADS-1) && gmesh->NElements()%NTHREADS) {
                todo[threadnumber]->push(gmesh->ElementVec()[(nelementsperthread*(threadnumber+1))]->Id());
            }
        } 
    }
	
	// end of initialization
	
	double time_start = get_time();
	cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS (V3)***\n";
    
	for (threadnumber=0; threadnumber<NTHREADS; threadnumber++)
    {
        divide_data *sdata;
        sdata = new divide_data;
        sdata->threadnumber = threadnumber;
        sdata->mesh = gmesh;
        threads.push_back(std::thread(divide,sdata));
    }
    
    for (threadnumber=0; threadnumber<NTHREADS; threadnumber++)
    {
        threads[threadnumber].join();
    }
	
	cout << "\n****END****\n";
	double time_end = get_time();
	cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
	
	ofstream af("after_PTHREAD3.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
	
	//bef.close();
	af.close();
    
	return 0;
}

#endif

//------------------------------------------------------------//

#ifdef REFTHREADSERIAL

TPZVec <TPZGeoEl *> sons;
double err;

int main ()
{ 	
	TPZGeoMesh *gmesh = GetMesh(NX, NY);
	
	double time_start = get_time();
	cout << "\n\n***STARTING SERIAL REFINEMENT PROCESS***\n";

	int iref, nref = MAXLVL;
	for (iref=0; iref<nref; iref++)
	{
		int nelements = gmesh->NElements();
		
		for (int iel=0; iel<nelements; iel++)
		{
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			
			#ifdef PZDEBUG
			if(!gel)
			{
				DebugStop();
			}
			#endif
			
			if(!gel->HasSubElement())
			{
				gel->Divide(sons);
			}
		}
	}
	
	cout << "\n****END****\n";
	double time_end = get_time();
	cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
	
	ofstream af("after_SERIAL.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
    
    af.close();
	
    //RunTestsCompMesh(gmesh);
    
	return 0;
}

#endif

//------------------------------------------------------------//

#ifdef ADMCHUNKTESTS

#include "pzadmchunkthreadsafe.h"

#define SIZEBLOCK 10
#define NITERATIONS NTHREADS*SIZEBLOCK

std::vector<std::thread> threads;

#if TEST==1
TPZAdmChunkVectorThreadSafe <int> cpyvect[NITERATIONS];
#endif

#if TEST==6
int objids[NITERATIONS];
#endif

struct threaddt
{
	TPZAdmChunkVectorThreadSafe <int> *vectin;
	int initid;
	int endid;
	int *obj[SIZEBLOCK];
	int threadid;
};

#if TEST==1
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	// int initid = internaldata->initid;
	// int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=0; i<NITERATIONS; i++)
	{
		cpyvect[i] = *internalvect;		// copying data
	}
	
}
#endif

#if TEST==2
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	// int initid = internaldata->initid;
	// int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=0; i<NITERATIONS; i++)
	{
		int ind = internalvect->AllocateNewElement();
		internalvect->operator[](ind) = ind;
	}
	
}
#endif

#if TEST==3
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	// int initid = internaldata->initid;
	// int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=0; i<NITERATIONS; i++)
	{
		internalvect->Resize(internalvect->NElements()+1);
	}
	
}
#endif

#if TEST==4
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	int initid = internaldata->initid;
	int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=initid; i<endid; i++)
	{
		internalvect->SetFree(i);
	}
	
}
#endif

#if TEST==5
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	int initid = internaldata->initid;
	int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=initid; i<endid; i++)
	{
		internalvect->operator[](i) = i;
	}
	
}
#endif

#if TEST==6
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	// int initid = internaldata->initid;
	// int endid = internaldata->endid;
	int threadid = internaldata->threadid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	int *obj[SIZEBLOCK];
	
	for(int i=0; i<SIZEBLOCK; i++)
	{
		obj[i] = internaldata->obj[i];
		
		objids[(threadid*SIZEBLOCK)+i] = internalvect->FindObject(obj[i]);
	}
	
}
#endif

#if TEST==7
void *function (void *arg)
{
	threaddt *internaldata = (threaddt *) arg;
	// int initid = internaldata->initid;
	// int endid = internaldata->endid;
	TPZAdmChunkVectorThreadSafe <int> *internalvect = internaldata->vectin;
	
	for (int i=0; i<NITERATIONS; i++)
	{
		internalvect->CompactDataStructure(2);
	}
	
}
#endif

int main()
{	
	int threadnumber;
	double err = 0.;
	TPZAdmChunkVectorThreadSafe <int> *vect = new TPZAdmChunkVectorThreadSafe <int> (0);
	
	threaddt thread_data[NTHREADS];
	
	for (int i=0; i<NTHREADS; i++)
	{
		thread_data[i].vectin = vect;
		
#if TEST==1
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=0; j<SIZEBLOCK; j++)
		{
			vect->AllocateNewElement();
		}
#endif
		
#if TEST==4
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=thread_data[i].initid; j<=thread_data[i].endid; j++)
		{
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
#endif
		
#if TEST==5
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=0; j<SIZEBLOCK; j++)
		{
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
#endif
		
#if TEST==6
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=thread_data[i].initid; j<=thread_data[i].endid; j++)
		{
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
		
		thread_data[i].threadid = i;
		
		for (int j=0; j<SIZEBLOCK; j++)
		{
			thread_data[i].obj[j] = &(thread_data[i].vectin->operator[]((i*SIZEBLOCK)+j));
		}
#endif
		
#if TEST==7
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=0; j<SIZEBLOCK; j++)
		{
			vect->AllocateNewElement();
		}
#endif
	}
	
	for (threadnumber=0; threadnumber<NTHREADS;)
	{
    threads.push_back(std::thread(function,&thread_data[threadnumber]));
    threadnumber++;
	}
	
#ifdef TEST
	for  (threadnumber=0; threadnumber<NTHREADS;)
	{
    threads[threadnumber].join();
    threadnumber++;
	}
#endif
	
	return 0;
}

#endif

//------------------------------------------------------------//
