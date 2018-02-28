#include <iostream>
#include <algorithm>
/**
 * PZ Includes
 */
#include "tpzdohrsubstruct.h"
#include "tpzdohrmatrix.h"
#include "tpzdohrprecond.h"
#include "pzdohrstructmatrix.h"
#include "pzdohrstructmatrix.cpp"
#include "pzstepsolver.h"
#include "pzcompel.h"
#include "pzgeoelbc.h"
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "tpzdohrassembly.h"
#include "pzlog.h"
#include "tpzgensubstruct.h"
#include "tpzpairstructmatrix.h"
#include "pzviscoelastic.h"
#include "TPZVTKGeoMesh.h"
#include "pzgeotetrahedra.h"
#include "pzskylstrmatrix.h"
#include "tpzarc3d.h"
#include "tpzdohrmatrix.h"
#include "pzvtkmesh.h"
#include "pzskylmat.h"

#ifdef USING_TBB
#include <tbb/blocked_range.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/parallel_for.h>
#endif

/**
 * Command Line Arguments
 */
clarg::argString    predio_file("-f", "Mesh file.", "8andares02.txt");
clarg::argInt       plevel("-p", "plevel", 1);
clarg::argInt       nsub("-nsub", "number of substructs", 32);
clarg::argInt       nloop("-l", "Number of loop iterations of the Subst_Backward/Subst_Forward", 1);
clarg::argBool      usetbb("-tbb", "Use of parallel tbb version", false);
clarg::argBool      aff_tbb("-aff", "Use of affinity partitioner", false);
clarg::argBool      help("-h", "Show usage message.", false);

RunStatsTable dec_rst   ("-dec", "Decompose statistics raw data table");
RunStatsTable sub_rst   ("-sub", "Substitution Forward/Backward statistics raw data table");


vector<TPZSkylMatrix<STATE>* > *get_sky_matrices();

#ifdef USING_TBB
class tbb_cholesky {
public:
    
    vector<TPZSkylMatrix<REAL>* > *fTasks;
    
    void operator()(const tbb::blocked_range<size_t>& range) const {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            (*fTasks)[i]->Decompose_Cholesky();
        }
    } // operator()
};

class tbb_substitution {
public:
    
    vector<TPZSkylMatrix<REAL>* > *fTasks;
    
    void operator()(const tbb::blocked_range<size_t>& range) const {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            TPZFMatrix<REAL> f((*fTasks)[i]->Dim(),1,M_PI);
            (*fTasks)[i]->Subst_Forward(&f);
            (*fTasks)[i]->Subst_Backward(&f);
        }
    } // operator()
};
#endif

void usage(char *prog)
{
    printf("\nUsage: %s\n", prog);
    printf("Arguments: \n");
    
    clarg::arguments_descriptions(cout, "   ", "\n");    
}

bool wayToSort(TPZSkylMatrix<REAL>* i, TPZSkylMatrix<REAL>* j) { return i->MemoryFootprint() > j->MemoryFootprint(); }

int main(int argc, char **argv)
{
    // parse the arguments
    if (clarg::parse_arguments(argc, argv)) {
        cerr << "Error when parsing the arguments!" << endl;
        usage(argv[0]);
        return 1;
    }
    
    if (!predio_file.was_set() || help.get_value()) {
        usage(argv[0]);
        return 1;
    }
    
    vector<TPZSkylMatrix<STATE>* > * fTasks = get_sky_matrices();
    
    std::sort(fTasks->begin(), fTasks->end(), wayToSort);
//    for(int i=0; i<fTasks->size(); i++)
//    {    
//        cout << (*fTasks)[i]->MemoryFootprint() << endl;// = new TPZSkylMatrix<REAL>(*orig);
//    }
    
    
    TPZSkylMatrix<REAL> *orig = (*fTasks)[0];
     fTasks->clear();
    fTasks->resize(nsub.get_value());
    for(int i=0; i<fTasks->size(); i++)
    {    
     (*fTasks)[i] = new TPZSkylMatrix<REAL>(*orig);
     }
     
     cout << "The copied matrix has Memory Foot Print of : " << (orig->MemoryFootprint()/(1024*1024)) << " MBs" << endl;
    
#ifdef USING_TBB
    tbb::task_scheduler_init init;
#endif
    
    int nmatrices = nsub.get_value();
    if (!usetbb.get_value()) {
        // serial decompose cholesky
        dec_rst.start();
        cout << "TPZSkylMatrix::Decompose_Cholesky()" << endl;
        for (int i=0; i<nmatrices; i++) {
            (*fTasks)[i]->Decompose_Cholesky();
        }
        dec_rst.stop();
        // serial Subst_Backward/Subst_Forward
        cout << "TPZSkylMatrix::Subst_Backward/Subst_Forward()" << endl;
        sub_rst.start();
        for (int k=0; k<nloop.get_value();k++) {
            for (int i=0; i<nmatrices; i++) {
                TPZFMatrix<REAL> f((*fTasks)[i]->Dim(),1,M_PI);
                (*fTasks)[i]->Subst_Forward(&f);
                (*fTasks)[i]->Subst_Backward(&f);
            }
        }
        sub_rst.stop();
    } else {
#ifdef USING_TBB
        
        tbb_cholesky disp;
        disp.fTasks = fTasks;
        
        tbb::affinity_partitioner ap;
         cout << "TPZSkylMatrix::Decompose_Cholesky()" << endl;
        dec_rst.start();
        if (aff_tbb.get_value())
            parallel_for(tbb::blocked_range<size_t>(0, nmatrices), disp, ap);
        else
            parallel_for(tbb::blocked_range<size_t>(0, nmatrices), disp);
        dec_rst.stop();
        tbb_substitution dispb;
        dispb.fTasks = fTasks;
        cout << "TPZSkylMatrix::Subst_Backward/Subst_Forward()" << endl;
        sub_rst.start();
        for (int k=0; k<nloop.get_value();k++) {
            if (aff_tbb.get_value())
                parallel_for(tbb::blocked_range<size_t>(0, nmatrices), dispb, ap);
            else
                parallel_for(tbb::blocked_range<size_t>(0, nmatrices), dispb);
        }
        sub_rst.stop();
#else
        cout << "Compiled without TBB support." << endl;
#endif
    }
    
    
    
} // main

TPZGeoMesh *malha_predio(string filename)
{
	//int nBCs = 1;
	int numnodes=-1;
	int numelements=-1;
	
	{
		bool countnodes = false;
		bool countelements = false;
		
		ifstream read (filename.c_str());
		if (!read.is_open()) {
            cerr << "Could not open file: " << filename << endl;
            exit(1);
		}
        
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "Coordinates") countnodes = true;
			if(str == "end coordinates") countnodes = false;
			if(countnodes) numnodes++;
			
			if(str == "Elements") countelements = true;
			if(str == "end elements") countelements = false;
			if(countelements) numelements++;
		}
	}
	
	TPZGeoMesh * gMesh = new TPZGeoMesh;
	
	gMesh -> NodeVec().Resize(numnodes);
	
	TPZVec <int64_t> TopolTetra(4);
	
	const int Qnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Qnodes);
	
	//setting nodes coords
	int64_t nodeId = 0, elementId = 0, matElId = 1;
	
	ifstream read;
	read.open(filename.c_str());
	
	double nodecoordX , nodecoordY , nodecoordZ ;
	
	char buf[1024];
	read.getline(buf, 1024);
	read.getline(buf, 1024);
	std::string str(buf);
	int in;
	for(in=0; in<numnodes; in++)
	{
		read >> nodeId;
		read >> nodecoordX;
		read >> nodecoordY;
		read >> nodecoordZ;
		Node[nodeId-1].SetNodeId(nodeId);
		Node[nodeId-1].SetCoord(0,nodecoordX);
		Node[nodeId-1].SetCoord(1,nodecoordY);
		Node[nodeId-1].SetCoord(2,nodecoordZ);
		gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
	}
	
	{
		read.close();
		read.open(filename.c_str());
        
		int l , m = numnodes+5;
		for(l=0; l<m; l++)
		{
			read.getline(buf, 1024);
		}
		
		int el;
		int matBCid = -1;
		//std::set<int> ncoordz; //jeitoCaju
		for(el=0; el<numelements; el++)
		{
			read >> elementId;
			read >> TopolTetra[0]; //node 1
			read >> TopolTetra[1]; //node 2
			read >> TopolTetra[2]; //node 3
			read >> TopolTetra[3]; //node 4
			
			// O GID comeca com 1 na contagem dos nodes, e nao zero como no PZ, assim o node 1 na verdade é o node 0
			TopolTetra[0]--;
			TopolTetra[1]--;
			TopolTetra[2]--;
			TopolTetra[3]--;
			
			int index = el;
			
			new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (index, TopolTetra, matElId, *gMesh);
		}
		
		gMesh->BuildConnectivity();
		// Colocando as condicoes de contorno
		
		for(el=0; el<numelements; el++)
		{
			TPZManVector <TPZGeoNode,4> Nodefinder(4);
			TPZManVector <REAL,3> nodecoord(3);
			TPZGeoEl *tetra = gMesh->ElementVec()[el];
			// na face z = 0
			TPZVec<int64_t> ncoordzVec(0); int64_t sizeOfVec = 0;
			for (int i = 0; i < 4; i++)
			{
				int64_t pos = tetra->NodeIndex(i);
				Nodefinder[i] = gMesh->NodeVec()[pos];
				Nodefinder[i].GetCoordinates(nodecoord);
                if (nodecoord[2] == 0.)
				{
					sizeOfVec++;
					ncoordzVec.Resize(sizeOfVec);
					ncoordzVec[sizeOfVec-1] = pos;
				}
			}
			if(ncoordzVec.NElements() == 3)
			{
				int lado = tetra->WhichSide(ncoordzVec);
				TPZGeoElSide tetraSide(tetra, lado);
				TPZGeoElBC(tetraSide,matBCid);
			}
		}
	}
	
	return gMesh;
}

void insert_elasticity(TPZAutoPointer<TPZCompMesh> mesh)
{
	mesh->SetDimModel(3);
	int nummat = 1;
	STATE E = 1.e6;
	STATE poisson = 0.3;
	TPZManVector<STATE> force(3,0.);
	force[1] = 20.;
	TPZElasticity3D *elast = new TPZElasticity3D(nummat,E,poisson,force);
	TPZMaterial * elastauto(elast);
	TPZFMatrix<STATE> val1(3,3,0.),val2(3,1,0.);
	TPZBndCond *bc = elast->CreateBC(elastauto, -1, 0, val1, val2);
	TPZMaterial * bcauto(bc);
	mesh->InsertMaterialObject(elastauto);
	mesh->InsertMaterialObject(bcauto);
}

template<class TVar>
class par_assemble_task_t
{
private:
    
    /** We divide the assemble procedure into N work items, which will
     be executed by one or several tasks. The TBB parallel_for
     construct automatically divide the work items in subsets and
     "creates" tasks to execute the work in each subset. Each task
     invokes the operator(blocked_range subset), which will be
     responsible for executing the work items in the subset. */
    template<class TTVar>
    struct work_item_t
    {
        work_item_t (unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TTVar> >& substruct) :
        fSubMeshIndex(submesh_idx), fSubstruct(substruct) {}
        
        unsigned fSubMeshIndex;
        TPZAutoPointer<TPZDohrSubstructCondense<TTVar> > fSubstruct;
    };
    
    /** Array of work items. */
    std::vector<work_item_t<TVar> > work_items;
    // TODO: Try the cache_aligned_allocator for improved performance.
    //std::vector<work_item_t<TVar>,cache_aligned_allocator<work_item_t<TVar> > > work_items;
    
    /* Pointers to shared data structures. */
    TPZAutoPointer<TPZDohrAssembly<TVar> > fAssembly;
    TPZAutoPointer<TPZCompMesh> fMesh;
    
public:
    
    par_assemble_task_t(TPZAutoPointer<TPZDohrAssembly<TVar> > assembly,
                        TPZAutoPointer<TPZCompMesh> mesh) :
    fAssembly(assembly), fMesh(mesh) {}
    
    /** Add a new work item to be list. */
    void push_work_item(unsigned submesh_idx, const TPZAutoPointer<TPZDohrSubstructCondense<TVar> >& substruct)
    {
        work_items.push_back(work_item_t<TVar>(submesh_idx, substruct));
    }
    
    /** Execute work items serially. */
    void run_serial()
    {
        typename std::vector<work_item_t<TVar> >::iterator it = work_items.begin();
        typename std::vector<work_item_t<TVar> >::iterator end = work_items.end();
        
        for (;it != end; it++)
        {
            work_item_t<TVar>& wi = *it;
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
        }
    }
    
#ifdef USING_TBB
    /** Computing operator for the parallel for. */
    void operator()(const blocked_range<size_t>& range) const
    {
        for(size_t i=range.begin(); i!=range.end(); ++i ) {
            const work_item_t<TVar>& wi = work_items[i];
            TPZSubCompMesh* submesh = SubMesh(fMesh, wi.fSubMeshIndex);
            AssembleMatrices(submesh, wi.fSubstruct, fAssembly,NULL);
        }
    }
    
    /** Execute work items in parallel. */
    void run_parallel_for()
    {
        /* TBB Parallel for. It will split the range into N sub-ranges and
         invoke the operator() for each sub-range.*/
        parallel_for(blocked_range<size_t>(0,work_items.size(), 1 /*IdealGrainSize*/), *this);
    }
#endif
};

vector<TPZSkylMatrix<STATE>* > * only_assemble(TPZDohrStructMatrix* dohrstruct,
                                               TPZMatrix<STATE> & mat, TPZFMatrix<STATE> & rhs,
                                               TPZAutoPointer<TPZCompMesh> &fMesh)
{
    TPZMatrix<STATE> *dohrgeneric = &mat;
    TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *dohr =
    dynamic_cast<TPZDohrMatrix<STATE,TPZDohrSubstructCondense<STATE> > *> (dohrgeneric);
    
    const std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > > &sublist = dohr->SubStructures();
    vector<TPZSkylMatrix<STATE>* > *fSkyTasks = new vector<TPZSkylMatrix<STATE>* >();
    std::list<TPZAutoPointer<TPZDohrSubstructCondense<STATE> > >::const_iterator it = sublist.begin();
    unsigned nsub = NSubMesh(fMesh);
    TPZSkylMatrix<STATE> *sky = 0;
    
    par_assemble_task_t<STATE> parallel_tasks(dohrstruct->Assembly(), fMesh);
    
    /* Initialize work items. */
    std::cout << "TPZDohrStructMatrix::Assemble()" << std::endl;
    for (unsigned isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
        if(!submesh) continue;
        parallel_tasks.push_work_item(isub, *it);
        it++;
    }
    
#ifdef USING_TBB
    parallel_tasks.run_parallel_for();
#else
    parallel_tasks.run_serial();
#endif
    
    it = sublist.begin();
    for (unsigned isub=0; isub<nsub ; isub++) {
        TPZSubCompMesh *submesh = SubMesh(fMesh, isub);
        if(!submesh) continue;
        
        // The Matred Big Matrix
        TPZAutoPointer<TPZMatRed<STATE,TPZFMatrix<STATE> > > matredbig = (*it)->fMatRedComplete;
        TPZAutoPointer<TPZMatrix<STATE> > Stiffness = matredbig->K00();
        
        sky = dynamic_cast<TPZSkylMatrix<STATE> *> (Stiffness.operator->());
        
        fSkyTasks->push_back(sky);
        
        // The Internal Matrix
        TPZAutoPointer<TPZMatRed<STATE,TPZFMatrix<STATE> > > matred = (*it)->fMatRed;
        TPZAutoPointer<TPZMatrix<STATE> > InternalStiffness = matred->K00();
        
        sky = dynamic_cast<TPZSkylMatrix<STATE> *> (InternalStiffness.operator->());
        
        fSkyTasks->push_back(sky);
        
        it++;
    }
    return fSkyTasks;
}

vector<TPZSkylMatrix<STATE>* > *get_sky_matrices() {
    
    TPZCompEl::SetgOrder(plevel.get_value());
    TPZGeoMesh  *gmesh = malha_predio(predio_file.get_value());
    TPZAutoPointer<TPZCompMesh> cmeshauto = new TPZCompMesh(gmesh);
    cmeshauto->SetDimModel(3);
    insert_elasticity(cmeshauto);
    cmeshauto->AutoBuild();
    
    TPZDohrStructMatrix* dohrstruct = new TPZDohrStructMatrix(cmeshauto);
    dohrstruct->IdentifyExternalConnectIndexes();
    
    cout << "TPZDohrStructMatrix::SubStructure()" << endl;
    dohrstruct->SubStructure(nsub.get_value());
    
    cout << "TPZDohrStructMatrix::Create()" << endl;
    TPZMatrix<STATE> *matptr = dohrstruct->Create();
    TPZFMatrix<STATE> *rhs = new TPZFMatrix<STATE>(cmeshauto->NEquations(),1,0.);
    
    return only_assemble(dohrstruct, *matptr,*rhs, cmeshauto);
}
