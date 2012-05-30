/// ---- GLOBALS ----

#define MAXLVL 5
#define NTHREADS 2

/// -----------

//#define REFTHREADTESTS			//enabling the tests
#define REFTHREAD					//main refthread
//#define REFTHREAD2					//main refthread2
//#define REFTHREADSERIAL			//serial implementation

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
#include "pzstack.h"
#include <pzgmesh.h>
#include <pzgeoel.h>
#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include <cmath>
#include "TPZRefPatternDataBase.h"

using namespace std;

pthread_t threads[NTHREADS];
pthread_mutex_t idslock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t todolock = PTHREAD_MUTEX_INITIALIZER;

TPZVec <int> *ids;
TPZStack <int> *todo;

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

TPZStack <int> getneighbours (TPZGeoMesh* mesh, int iel) 
{
	TPZStack <int> nbs;
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
				nneighbours = nbs.size();
				for (int j=0; j<nneighbours; j++) {
					if (neighbour.Element()->Id()!=nbs[j]) countneighbours++;
				}
				if (countneighbours==nneighbours) {
					int pos = neighbour.Element()->Id();
					nbs.Push(pos);
				}
				neighbour = neighbour.Neighbour();
			}
		}
	}
	return nbs;
}

bool checkneighbours (TPZGeoMesh* mesh, int iel) 
{
	if (ids->operator[](iel)) {
		return true;
	}
	
	else {
		int countneighbours=0;
		int nneighbours=0;
		
		TPZStack <int> neighbours = getneighbours(mesh, iel);
		
		nneighbours = neighbours.size();
		
		while (neighbours.size()) {
			int neigh = neighbours.Pop();
			if(!ids->operator[](neigh)) countneighbours++;
		}
		
		if (countneighbours==nneighbours) {
			return false;
		}
		
		else {
			return true;
		}
	}
}

int lockneighbours (TPZGeoMesh *mesh, int iel)
{
	pthread_mutex_lock(&idslock);
	
	if (checkneighbours(mesh, iel)) {
		pthread_mutex_lock(&todolock);
		todo->Push(iel);
		pthread_mutex_unlock(&todolock);
	
		pthread_mutex_unlock(&idslock);
		
		pthread_exit(0);
	}
	
	else {
		TPZStack <int> neighbours = getneighbours(mesh, iel);
		
		while (neighbours.size()) {
			int pos = neighbours.Pop();
			ids->Fill((int) 1, pos, 1);
		}
		
		ids->Fill((int) 1, iel, 1);
		
		pthread_mutex_unlock(&idslock);
		return 0;
	}
}

int unlockneighbours (TPZGeoMesh *mesh, int iel, TPZVec <TPZGeoEl*> sons)
{
	pthread_mutex_lock(&idslock);

	TPZStack <int> neighbours = getneighbours(mesh, iel);
	
	ids->Fill(0, iel, 1);
	
	while (neighbours.size()) {
		int pos = neighbours.Pop();
		ids->Fill(0, pos, 1);
	}
	
	pthread_mutex_unlock(&idslock);
	
	pthread_mutex_lock(&todolock);
	for (int i=0; i<sons.NElements(); i++) {
		if (sons[i]->Level()<MAXLVL) {
			todo->Push(sons[i]->Id());
		}
	}
	pthread_mutex_unlock(&todolock);
	
	return 0;
}

void *divide(void *arg)
{
	divide_data *idata = (divide_data*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	int iel = idata->el;
	delete idata;
	
	TPZVec <TPZGeoEl *> sons;
	
	TPZGeoEl * gel = imesh->ElementVec()[iel];
	
#ifdef DEBUG
	if(!gel)
	{
		DebugStop();
	}
#endif
	
	lockneighbours(imesh, iel);
	
	gel->Divide(sons);
	
	unlockneighbours(imesh, iel, sons);
	
	pthread_exit(0);
}

int main ()
{ 
	pthread_mutex_init(&idslock, NULL);
	pthread_mutex_init(&todolock, NULL);
	
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	
	REAL coordinates [16][3] = {{0.,0.,0.},{1.,0.,0.},{2.,0.,0.},{3.,0.,0.},
								{0.,1.,0.},{1.,1.,0.},{2.,1.,0.},{3.,1.,0.},
								{0.,2.,0.},{1.,2.,0.},{2.,2.,0.},{3.,2.,0.},
								{0.,3.,0.},{1.,3.,0.},{2.,3.,0.},{3.,3.,0.}};	// Instancing the coordinates
	int nodeids [16] = {100,1100,200,300,
						220,111,223,111,
						333,444,555,666,
						111,222,333,44};	// Identifying the nodes
	int elconnect[9][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{4,5,9,8},{5,6,10,9},{6,7,11,10},{8,9,13,12},{9,10,14,13},{10,11,15,14}};
	int elid [9] = {10,11,12,13,14,15};
	
	int index = 0;
	int i, j;	// Auxiliary
	
	TPZVec<REAL> coord(3,0.);	// Creating a vector of REAL (based on template) with size 3, initialized with zeros (0.)
	
	for (i=0; i<16; i++)
	{
		for (j=0; j<3; j++) coord[j] = coordinates[i][j];		// Going through the coordinates
		
		index = gmesh->NodeVec().AllocateNewElement();		// Allocating space for new elements (returned the index of the next free element or of the new element created)
		gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);		// Constructing the nodes, with the id got from the vector nodeids, with the TPZVec coord, on the mesh gmesh;
	}
	
	TPZGeoEl *elvec[9];		// Creating a vector of TPZGeoEl using the default constructor;
	
	TPZVec<int> connect(4,0);	// Auxiliary TPZVec that will be used on CreateGeoElement;
	
	for (i=0; i<9; i++)
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
	
	int threadnumber = 0;
	int totalElements = 1;
	REAL alevel;
	
	for (alevel=1; alevel<=MAXLVL; alevel++) {
		totalElements+=pow(4., alevel);
	}
	
	totalElements*=gmesh->NElements();
	
	ids = new TPZVec <int> (totalElements, azero);
	
	todo = new TPZStack <int>;
	
	for (int i=0; i<gmesh->NElements(); i++) {
		todo->Push(gmesh->ElementVec()[i]->Id());
	}
	
	while (todo->size()||threadnumber)
	{
		int iel = -1;
		
		if (todo->size()) {
			pthread_mutex_lock(&todolock);
			iel = todo->Pop();
			pthread_mutex_unlock(&todolock);
		}
		
		divide_data *sdata;
		sdata = new divide_data;
		sdata->el = iel;
		sdata->mesh = gmesh;
		
		if (threadnumber==NTHREADS||!todo->size())
		{
			int limit=threadnumber;
			
			int i = 0;
			for (;i<limit;)
			{
				err = pthread_join (threads[i], NULL);
				if (err)
				{
					cout << "Could not join the threads! Return code from pthread_join is " << err;
					DebugStop();
				}
				else
				{
					i++;
				}
			}
			
			threadnumber = 0;
		}
		if (iel!=-1) { 
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

#ifdef REFTHREAD2

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

pthread_t threads[NTHREADS];

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

void *divide(void *arg)
{
	divide_data *idata = (divide_data*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	int iel = idata->el;
	delete idata;
	
	TPZVec <TPZGeoEl *> sons;
	
	TPZStack <TPZGeoEl *> todo;
	
	todo.Push(imesh->ElementVec()[iel]);
	
	while (todo.size()) {
		TPZGeoEl * gel = todo.Pop();
		
		gel->Divide(sons);
	
		for (int i=0; i<sons.NElements(); i++) {
			if (sons[i]->Level()<MAXLVL) {
				todo.Push(sons[i]);
			}
		}
	}
	
	pthread_exit(0);
}

int main ()
{ 
	TPZGeoMesh *gmesh = new TPZGeoMesh;
	
	REAL coordinates [16][3] = {{0.,0.,0.},{1.,0.,0.},{2.,0.,0.},{3.,0.,0.},
		{0.,1.,0.},{1.,1.,0.},{2.,1.,0.},{3.,1.,0.},
		{0.,2.,0.},{1.,2.,0.},{2.,2.,0.},{3.,2.,0.},
		{0.,3.,0.},{1.,3.,0.},{2.,3.,0.},{3.,3.,0.}};	// Instancing the coordinates
	int nodeids [16] = {100,1100,200,300,
		220,111,223,111,
		333,444,555,666,
		111,222,333,44};	// Identifying the nodes
	int elconnect[9][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{4,5,9,8},{5,6,10,9},{6,7,11,10},{8,9,13,12},{9,10,14,13},{10,11,15,14}};
	int elid [9] = {10,11,12,13,14,15};
	
	int index = 0;
	int i, j;	// Auxiliary
	
	TPZVec<REAL> coord(3,0.);	// Creating a vector of REAL (based on template) with size 3, initialized with zeros (0.)
	
	for (i=0; i<16; i++)
	{
		for (j=0; j<3; j++) coord[j] = coordinates[i][j];		// Going through the coordinates
		
		index = gmesh->NodeVec().AllocateNewElement();		// Allocating space for new elements (returned the index of the next free element or of the new element created)
		gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);		// Constructing the nodes, with the id got from the vector nodeids, with the TPZVec coord, on the mesh gmesh;
	}
	
	TPZGeoEl *elvec[9];		// Creating a vector of TPZGeoEl using the default constructor;
	
	TPZVec<int> connect(4,0);	// Auxiliary TPZVec that will be used on CreateGeoElement;
	
	for (i=0; i<9; i++)
	{
		for (j=0; j<4; j++) connect[j] = elconnect[i][j];
		elvec[i] = gmesh->CreateGeoElement(EQuadrilateral, connect, 1, elid[i]);		// Creating the elements with type EQuadrilateral, node indexes from connect, with the material with id=1, indexes from vector elid
		// Returns on the vector of TPZGeoEl elvec
	}
	
	gmesh->BuildConnectivity();		// Building connectivity on the mesh
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
	
	/*for (int iref=0; iref<2; iref++) //initializing mesh for
	{
		int nelements = gmesh->NElements();
		TPZVec <TPZGeoEl *> sons;
		
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
	}*/
	
	ofstream bef("before_PTHREAD-2.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file
	
	double time_start = get_time();
	cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS***\n";
	
	int threadnumber = 0;
	
	int todoelements[4] = {0, 2, 6, 8};

	for  (threadnumber = 0; threadnumber<NTHREADS;)
	{
		divide_data *sdata;
		sdata = new divide_data;
		sdata->el = todoelements[threadnumber];
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
	
	for  (threadnumber = 0; threadnumber<NTHREADS;)
	{
		err = pthread_join (threads[threadnumber], NULL);
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
	
	cout << "\n****END****\n";
	double time_end = get_time();
	cout << "This refinement has run for: "<< (time_end-time_start) << " seconds.\n\n";
	
	ofstream af("after_PTHREAD-2.vtk");
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
	
	REAL coordinates [16][3] = {{0.,0.,0.},{1.,0.,0.},{2.,0.,0.},{3.,0.,0.},
		{0.,1.,0.},{1.,1.,0.},{2.,1.,0.},{3.,1.,0.},
		{0.,2.,0.},{1.,2.,0.},{2.,2.,0.},{3.,2.,0.},
		{0.,3.,0.},{1.,3.,0.},{2.,3.,0.},{3.,3.,0.}};	// Instancing the coordinates
	int nodeids [16] = {100,1100,200,300,
		220,111,223,111,
		333,444,555,666,
		111,222,333,44};	// Identifying the nodes
	int elconnect[9][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{4,5,9,8},{5,6,10,9},{6,7,11,10},{8,9,13,12},{9,10,14,13},{10,11,15,14}};
	int elid [9] = {10,11,12,13,14,15};
	
	int index = 0;
	int i, j;	// Auxiliary
	
	TPZVec<REAL> coord(3,0.);	// Creating a vector of REAL (based on template) with size 3, initialized with zeros (0.)
	
	for (i=0; i<16; i++)
	{
		for (j=0; j<3; j++) coord[j] = coordinates[i][j];		// Going through the coordinates
		
		index = gmesh->NodeVec().AllocateNewElement();		// Allocating space for new elements (returned the index of the next free element or of the new element created)
		gmesh->NodeVec()[index] = TPZGeoNode(nodeids[i],coord,*gmesh);		// Constructing the nodes, with the id got from the vector nodeids, with the TPZVec coord, on the mesh gmesh;
	}
	
	TPZGeoEl *elvec[9];		// Creating a vector of TPZGeoEl using the default constructor;
	
	TPZVec<int> connect(4,0);	// Auxiliary TPZVec that will be used on CreateGeoElement;
	
	for (i=0; i<9; i++)
	{
		for (j=0; j<4; j++) connect[j] = elconnect[i][j];
		elvec[i] = gmesh->CreateGeoElement(EQuadrilateral, connect, 1, elid[i]);		// Creating the elements with type EQuadrilateral, node indexes from connect, with the material with id=1, indexes from vector elid
		// Returns on the vector of TPZGeoEl elvec
	}
	
	gmesh->BuildConnectivity();		// Building connectivity on the mesh
	gRefDBase.InitializeUniformRefPattern(EQuadrilateral);	
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
	
	ofstream af("after_SERIAL.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
	
	return 0;
}

#endif

