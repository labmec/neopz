/// -----------

//#define REFTHREADTESTS			//enabling the tests
#define REFTHREAD				//main refthread
//#define REFTHREADSERIAL		//serial implementation

/// On REFTHREADTESTS:

#define TEST 7

// if TEST = 1				//for testing operator=() on the vector
// if TEST = 2				//for testing AllocateNewElement() on the vector
// if TEST = 3				//for testing Resize() on the vector
// if TEST = 4				//for testing SetFree() on the vector
// if TEST = 5				//for testing operator[]() on the vector
// if TEST = 6				//for testing FindObject() on the vector
// if TEST = 7				//for testing CompactDataStructure() on the vector

/// -----------

#ifdef REFTHREADTESTS

#include "pzadmchunkthreadsafe.h"

#include <iostream>
#include <pthread.h>
#include <sys/time.h>

#define NTHREADS 100
#define SIZEBLOCK 10
#define NITERATIONS NTHREADS*SIZEBLOCK

using namespace std;

pthread_t threads [NTHREADS];

#if TEST==1
TPZAdmChunkVectorThreadSafe <int> cpyvect[NITERATIONS];
#endif

#if TEST==6
int objids[NITERATIONS];
#endif

struct threaddt {
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
	
	pthread_exit(0);
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
	
	pthread_exit(0);
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
	
	pthread_exit(0);
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
	
	pthread_exit(0);
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
	
	pthread_exit(0);
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
	
	for(int i=0; i<SIZEBLOCK; i++) {
		obj[i] = internaldata->obj[i];
		
		objids[(threadid*SIZEBLOCK)+i] = internalvect->FindObject(obj[i]);
	}
	
	pthread_exit(0);
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
	
	pthread_exit(0);
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
		
		for (int j=0; j<SIZEBLOCK; j++) {
			vect->AllocateNewElement();
		}
		#endif
		
		#if TEST==4
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=thread_data[i].initid; j<=thread_data[i].endid; j++) {
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
		#endif
			
		#if TEST==5
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=0; j<SIZEBLOCK; j++) {
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
		#endif
		
		#if TEST==6
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=thread_data[i].initid; j<=thread_data[i].endid; j++) {
			vect->AllocateNewElement();
			vect->operator[](j) = j;
		}
		
		thread_data[i].threadid = i;
		
		for (int j=0; j<SIZEBLOCK; j++) {
			thread_data[i].obj[j] = &(thread_data[i].vectin->operator[]((i*SIZEBLOCK)+j));
		}
		#endif
		
		#if TEST==7
		thread_data[i].initid = (i*SIZEBLOCK);
		thread_data[i].endid = (i*SIZEBLOCK + (SIZEBLOCK-1));
		
		for (int j=0; j<SIZEBLOCK; j++) {
			vect->AllocateNewElement();
		}
		#endif
	}
	
	for (threadnumber=0; threadnumber<NTHREADS;)
	{
		err = pthread_create (&threads[threadnumber], NULL, function, (void *) &thread_data[threadnumber]);
		
		if (err)
		{
			cout << "There is a problem on creating the thread! Exiting the program! Return code from pthread_create is " << err;
			DebugStop();
		}
		else
		{
			threadnumber++;
		}
	}
	
	#ifdef TEST
	for  (threadnumber=0; threadnumber<NTHREADS;)
	{
		err = pthread_join (threads[threadnumber], NULL);
		if (err)
		{
			cout << "Could not join the threads! Return code from pthread_join is " << err;
			DebugStop();
		}
		else
		{
			threadnumber++;
		}
	}
	#endif
	
	return 0;
}

#endif

//------------------------------------------------------------//

#ifdef REFTHREAD

#include "pzvec.h"
#include <pzgmesh.h>
#include <pzgeoel.h>
#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include <cmath>
#include "TPZRefPatternDataBase.h"

using namespace std;

#define MAXLVL 10
#define NTHREADS 6
#define INTTHREAD 0

pthread_t threads[NTHREADS];
pthread_mutex_t datalock = PTHREAD_MUTEX_INITIALIZER;

TPZVec <int> *ids;

double err;
const int azero = 0;

struct divide_data {
	int el;
	TPZGeoMesh* mesh;
};

double get_time()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	double d = t.tv_sec + (double) t.tv_usec/1000000;
	return d;
}	// Simple function for taking time

TPZVec <int> getneighbours (TPZGeoMesh* mesh, int iel) 
{
	TPZVec <int> nbs;
	int countneighbours=0, nneighbours=0;
	
	int nsides = mesh->ElementVec()[iel]->NSides();
	for (int i=0; i<nsides; i++) {
		TPZGeoElSide neighbour = mesh->ElementVec()[iel]->Neighbour(i);
		TPZGeoElSide thisside (mesh->ElementVec()[iel],i);
		if (neighbour.Exists()) 
		{
			int count=0;
			while (neighbour!=thisside&&count++<30) {
				countneighbours = 0;
				nneighbours = nbs.NElements();
				for (int j=0; j<nneighbours; j++) {
					if (neighbour.Element()->Id()!=nbs[j]) countneighbours++;
				}
				if (countneighbours==nneighbours) nbs.Resize(nbs.NElements()+1, neighbour.Element()->Id());
				neighbour = neighbour.Neighbour();
			}
		}
	}
	return nbs;
}

int checkneighbours (TPZGeoMesh* mesh, int iel) 
{
	int countneighbours=0;
	int nneighbours=0;
	
	TPZVec <int> neighbours = getneighbours(mesh, iel);

	nneighbours = neighbours.NElements();
	
	for (int j=0; j<nneighbours; j++) {
		pthread_mutex_lock(&datalock);
		int pos = neighbours[j];
		int value = ids->operator[](pos);
		if(value==0) countneighbours++;
		pthread_mutex_unlock(&datalock);
		
	}
	
	if (countneighbours==nneighbours) return 0;
	else return 1;
}

int lockneighbours (TPZGeoMesh *mesh, int iel)
{
	int nneighbours=0;
	
	TPZVec <int> neighbours = getneighbours(mesh, iel);
	
	nneighbours = neighbours.NElements();
	for (int j=0; j<nneighbours; j++) {
		pthread_mutex_lock(&datalock);
		ids->Fill(1, neighbours[j], 1);
		pthread_mutex_unlock(&datalock);
	}
	pthread_mutex_lock(&datalock);
	ids->Fill(1, iel, 1);
	pthread_mutex_unlock(&datalock);
	
	return 0;
}

int unlockneighbours (TPZGeoMesh *mesh, int iel)
{
	int nneighbours=0;
	
	TPZVec <int> neighbours = getneighbours(mesh, iel);

	nneighbours = neighbours.NElements();
	
	for (int j=0; j<nneighbours; j++) {
		pthread_mutex_lock(&datalock);
		ids->Fill(0, neighbours[j], 1);
		pthread_mutex_unlock(&datalock);
	}
	pthread_mutex_lock(&datalock);
	ids->Fill(0, iel, 1);
	pthread_mutex_unlock(&datalock);
	
	return 0;
}

void *divide(void *arg)
{
	divide_data *idata = (divide_data*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	int iel = idata->el;
	delete idata;
	//cout << "---Refining element no: " <<  iel << ". ----\n";
	
	TPZVec <TPZGeoEl *> sons;
	
	TPZGeoEl * gel = imesh->ElementVec()[iel];
	
#ifdef DEBUG
	if(!gel)
	{
		DebugStop();
	}
#endif
	
	err = lockneighbours(imesh, iel);
	
	if(!err) gel->Divide(sons);
	
	err = unlockneighbours(imesh, iel);
	
	if(err) {
		DebugStop();
	}
	
	pthread_exit(0);
}

int main ()
{ 
	pthread_mutex_init(&datalock, NULL);
	
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	
	REAL coordinates [6][3] = {{0.,0.,0.},{1.,0.,1.},{2.,0.,0.},{0.,1.,0.},{1.,1.,1.},{2.,1.,0.}};	// Instancing the coordinates
	int nodeids [6] = {100,1100,200,300,1200,400};	// Identifying the nodes
	int elconnect[2][4] = {{0,1,4,3},{4,1,2,5}};
	int elid [2] = {10,20};
	
	int index = 0;
	int i, j;	// Auxiliary
	
	TPZVec<REAL> coord(3,0.);	// Creating a vector of REAL (based on template) with size 3, initialized with zeros (0.)
	
	for (i=0; i<6; i++)
	{
		for (j=0; j<3; j++) coord[j] = coordinates[i][j];		// Going through the coordinates
		
		index = gmesh->NodeVec().AllocateNewElement();		// Allocating space for new elements (returned the index of the next free element or of the new element created)
		gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);		// Constructing the nodes, with the id got from the vector nodeids, with the TPZVec coord, on the mesh gmesh;
	}
	
	TPZGeoEl *elvec[2];		// Creating a vector of TPZGeoEl using the default constructor;
	
	TPZVec<int> connect(4,0);	// Auxiliary TPZVec that will be used on CreateGeoElement;
	
	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++) connect[j] = elconnect[i][j];
		elvec[i] = gmesh->CreateGeoElement(EQuadrilateral, connect, 1, elid[i]);		// Creating the elements with type EQuadrilateral, node indexes from connect, with the material with id=1, indexes from vector elid
		// Returns on the vector of TPZGeoEl elvec
	}
	
	gmesh->BuildConnectivity();		// Building connectivity on the mesh
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	
	ofstream bef("before_PTHREAD.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file
	
	double time_start = get_time();
	cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS***\n";
	
	int threadnumber = INTTHREAD;
	int totalElements = 1;
	REAL alevel;
	
	for (alevel=1; alevel<=MAXLVL; alevel++) {
		totalElements+=pow(4., alevel);
	}
	
	totalElements*=gmesh->NElements();

	ids = new TPZVec <int> (totalElements, azero);
	
	do
	{
		int nelements = gmesh->NElements();
		for (int iel=0; iel<nelements; iel++)
		{	
			if (!ids->operator[](iel) && gmesh->ElementVec()[iel]->Level()<MAXLVL && !gmesh->ElementVec()[iel]->HasSubElement() && !checkneighbours(gmesh, iel))
			{
				if (threadnumber==NTHREADS)
				{
					threadnumber = INTTHREAD;
					for  (;threadnumber<NTHREADS;)
					{
						err = pthread_join (threads[threadnumber], NULL);
						if (err)
						{
							cout << "Could not join the threads! Return code from pthread_join is " << err;
							DebugStop();
						}
						else
						{
							threadnumber++;
						}
					}
					threadnumber = INTTHREAD;
				}
				
				divide_data *sdata;
				sdata = new divide_data;
				sdata->el = iel;
				sdata->mesh = gmesh;
				
				err = pthread_create (&threads[threadnumber], NULL, divide, (void *) sdata);
				if (err)
				{
					cout << "There is a problem on creating the thread! Exiting the program! Return code from pthread_create is " << err;
					DebugStop();
				}
				else
				{
					threadnumber++;
					
				}
			}
		}
	} while (gmesh->NElements()<totalElements);
	
	cout << "\n****END****\n";
	double time_end = get_time();
	cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
	
	ofstream af("after_PTHREAD.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
	
	bef.close();
	af.close();
			
	return 0;
}

#endif

//------------------------------------------------------------//

#ifdef REFTHREADSERIAL

#include "pzvec.h"
#include <pzgmesh.h>
#include <pzgeoel.h>
#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include "TPZRefPatternDataBase.h"

using namespace std;

TPZVec <TPZGeoEl *> sons;
double err;

#define MAXLVL 3

double get_time()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	double d = t.tv_sec + (double) t.tv_usec/1000000;
	return d;
}	// Simple function for taking time

int main ()
{ 	
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	
	REAL coordinates [6][3] = {{0.,0.,0.},{1.,0.,1.},{2.,0.,0.},{0.,1.,0.},{1.,1.,1.},{2.,1.,0.}};	// Instancing the coordinates
	
	int nodeids [6] = {100,1100,200,300,1200,400};	// Identifying the nodes
	
	int elconnect[2][4] = {{0,1,4,3},{4,1,2,5}};
	
	int elid [2] = {10,20};
	
	int index = 0;
	int i, j;	// Auxiliary

	TPZVec<REAL> coord(3,0.);	// Creating a vector of REAL (based on template) with size 3, initialized with zeros (0.)
	
	for (i=0; i<6; i++)
	{
		for (j=0; j<3; j++) coord[j] = coordinates[i][j];		// Going through the coordinates
		
		index = gmesh->NodeVec().AllocateNewElement();		// Allocating space for new elements (returned the index of the next free element or of the new element created)
		gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);		// Constructing the nodes, with the id got from the vector nodeids, with the TPZVec coord, on the mesh gmesh;
	}
	
	TPZGeoEl *elvec[2];		// Creating a vector of TPZGeoEl using the default constructor;
	
	TPZVec<int> connect(4,0);	// Auxiliary TPZVec that will be used on CreateGeoElement;
	
	for (i=0; i<2; i++)
	{
		for (j=0; j<4; j++) connect[j] = elconnect[i][j];
		elvec[i] = gmesh->CreateGeoElement(EQuadrilateral, connect, 1, elid[i]);		// Creating the elements with type EQuadrilateral, node indexes from connect, with the material with id=1, indexes from vector elid
		// Returns on the vector of TPZGeoEl elvec
	}

	gmesh->BuildConnectivity();		// Building connectivity on the mesh

	
	//gmesh->Print();					// Printing the initial mesh
	ofstream bef("before_SERIAL.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file

	double time_start = get_time();
	cout << "\n\n***STARTING SERIAL REFINEMENT PROCESS***\n";
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);

	int iref, nref = MAXLVL;
	for (iref=0; iref<nref; iref++)
	{
		int nelements = gmesh->NElements();
		
		for (int iel=0; iel<nelements; iel++)
		{
			TPZGeoEl * gel = gmesh->ElementVec()[iel];
			
			#ifdef DEBUG
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
	
	//gmesh->Print();					// Printing the initial mesh
	ofstream af("after_SERIAL.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
	
	return 0;
}

#endif

