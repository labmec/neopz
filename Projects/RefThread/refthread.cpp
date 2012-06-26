/// ---- GLOBALS ----

#define MAXLVL 5
#define NTHREADS 2
#define NX 5
#define NY 5

/// -----------

//#define REFTHREADV2					//main refthread second version
#define REFTHREADV1						//main refthread first working version
//#define REFTHREADSERIAL				//serial implementation
//#define REFTHREADTESTS				//main refthreadtests
//#define ADMCHUNKTESTS					//enabling the tests on admchunkvectorthreadsafe

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

#include <iostream>
#include <pthread.h>
#include <sys/time.h>
#include <cmath>

#include <pzgmesh.h>
#include <pzgeoel.h>
#include "TPZVTKGeoMesh.h"
#include "TPZRefPatternDataBase.h"

using namespace std;

double get_time()
{
	struct timeval t;
	gettimeofday(&t, NULL);
	double d = t.tv_sec + (double) t.tv_usec/1000000;
	return d;
}	// Simple function for taking time

TPZGeoMesh* GetMesh (int nx, int ny){
	int i,j;
	int id, index;
	
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
	
	TPZVec <int> connect(4,0);
	
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

//------------------------------------------------------------//

#ifdef REFTHREADV2

#include "pzstack.h"
#include <set>

pthread_t threads[NTHREADS];
pthread_mutex_t todolock = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t elsetlock = PTHREAD_MUTEX_INITIALIZER;

TPZStack <int> *todo;
set <int> elset;

double err;
const int azero = 0;

struct divide_data 
{
	int el;
	TPZGeoMesh* mesh;
};

TPZStack <int> getneighbours (TPZGeoMesh* mesh, int iel) 
{
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
	return nbs;
}

bool checkneighbours (TPZGeoMesh* mesh, int iel) 
{	
	int found = (*(elset.find(iel))==iel);
	if (!found)
	{
		return false;
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
			found = !(elset.find(neigh)==elset.end());
			if(found) countneighbours++;
		}
		
		if (countneighbours==nneighbours)
		{
			return true;
		}
		
		else
		{
			return false;
		}
	}
}

int lockneighbours (TPZGeoMesh *mesh, int iel)
{
	
	if (checkneighbours(mesh, iel)) {
		
		elset.insert(iel);
		
		pthread_mutex_unlock(&elsetlock);
		
		return 1;
	}
	
	else {
		TPZStack <int> neighbours = getneighbours(mesh, iel);
		
		while (neighbours.size())
		{
			int pos = neighbours.Pop();
			elset.erase(pos);
		}
		
		pthread_mutex_unlock(&elsetlock);
		
		return 0;
	}
}

int unlockneighbours (TPZGeoMesh *mesh, int iel, TPZVec <TPZGeoEl*> sons)
{
	pthread_mutex_lock(&elsetlock);
	
	TPZStack <int> neighbours = getneighbours(mesh, iel);
	
	elset.insert(iel);
	
	while (neighbours.size()) {
		int pos = neighbours.Pop();
		elset.insert(pos);
	}
	
	for (int i=0; i<sons.NElements(); i++) {
		if (sons[i]->Level()<MAXLVL) {
			elset.insert(sons[i]->Id());
		}
	}
	
	pthread_mutex_unlock(&elsetlock);
	return 0;
}

void *divide(void *arg)
{
	divide_data *idata = (divide_data*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	int iel = idata->el;
	delete idata;
	
	TPZVec <TPZGeoEl *> sons;
	
	imesh->ElementVec()[iel]->Divide(sons);
	
	unlockneighbours(imesh, iel, sons);
	
	pthread_exit(0);
}

int main ()
{ 
	pthread_mutex_init(&todolock, NULL);
	pthread_mutex_init(&elsetlock, NULL);
	
	TPZGeoMesh *gmesh = GetMesh(NX, NY);
	
	ofstream bef("before_PTHREAD.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file
	
	double time_start = get_time();
	cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS***\n";
	
	int threadnumber = 0;
	
	for (int i=0; i<gmesh->NElements(); i++) {
		elset.insert(gmesh->ElementVec()[i]->Id());
	}
	
	while (elset.size()||threadnumber)
	{
		int iel = -1;
		
		if (elset.size()) {
			pthread_mutex_lock(&elsetlock);
			set<int>::iterator it = elset.begin();
			advance(it, rand() % elset.size());
			iel = *(it);
			elset.erase(it);
			
			int err = lockneighbours (gmesh, iel);
			if (err)
			{
				continue;
			}
		}
		
		divide_data *sdata;
		sdata = new divide_data;
		sdata->el = iel;
		sdata->mesh = gmesh;
		
		if (threadnumber==NTHREADS||!elset.size())
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

#ifdef REFTHREADV1

#include "pzvec.h"
#include "pzstack.h"

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

TPZStack <int> getneighbours (TPZGeoMesh* mesh, int iel) 
{
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
	return nbs;
}

bool checkneighbours (TPZGeoMesh* mesh, int iel) 
{
	if (ids->operator[](iel))
	{
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
			return false;
		}
		
		else
		{
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
	
	TPZGeoMesh *gmesh = GetMesh(NX, NY);
	
	ofstream bef("before_PTHREAD.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, bef);			// Printing the initial mesh on the file
	
	double time_start = get_time();
	cout << "\n\n***STARTING PTHREAD REFINEMENT PROCESS***\n";
	
	int threadnumber = 0;
	REAL totalElements = 1;
	REAL alevel;
	
	for (alevel=1; alevel<=MAXLVL; alevel++)
	{
		totalElements+=pow(4., alevel);
	}
	
	totalElements*=gmesh->NElements();
	
	ids = new TPZVec <int> (totalElements, azero);
	
	todo = new TPZStack <int>;
	
	for (int i=0; i<gmesh->NElements(); i++)
	{
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
		
		if (iel!=-1)
		{
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

#ifdef REFTHREADSERIAL

#include "pzvec.h"

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

//------------------------------------------------------------//

#ifdef REFTHREADTESTS

#include "pzvec.h"

pthread_t threads[NTHREADS];

double err;
const int azero = 0;

struct divide_data
{
	int el;
	TPZGeoMesh* mesh;
};

void *divide(void *arg)
{
	divide_data *idata = (divide_data*) arg;
	TPZGeoMesh *imesh = idata->mesh;
	int iel = idata->el;
	delete idata;
	
	TPZVec <TPZGeoEl *> sons;
	
	TPZStack <TPZGeoEl *> todo;
	
	todo.Push(imesh->ElementVec()[iel]);
	
	while (todo.size())
	{
		TPZGeoEl * gel = todo.Pop();
		
		gel->Divide(sons);
		
		for (int i=0; i<sons.NElements(); i++)
		{
			if (sons[i]->Level()<MAXLVL)
			{
				todo.Push(sons[i]);
			}
		}
	}
	
	pthread_exit(0);
}

int main ()
{ 
	TPZGeoMesh *gmesh = GetMesh(NX, NY);
	
	ofstream bef("before_PTHREAD_tests.vtk");
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
	
	ofstream af("after_PTHREAD_tests.vtk");
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, af);
	
	bef.close();
	af.close();
	
	return 0;
}

#endif

//------------------------------------------------------------//

#ifdef ADMCHUNKTESTS

#include "pzadmchunkthreadsafe.h"

#define SIZEBLOCK 10
#define NITERATIONS NTHREADS*SIZEBLOCK

pthread_t threads [NTHREADS];

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
	
	for(int i=0; i<SIZEBLOCK; i++)
	{
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