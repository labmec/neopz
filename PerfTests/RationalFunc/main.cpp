/**
 * @file
 * @brief Projeto elaborado para avaliar o desempenho paralelo da matriz frontal.
 */


#include "arglib.h"             // clarg:: command line arguments lib
#include "run_stats_table.h"
#include "CedricTest.h"         // TCedricTest
#include "TPZRefPatternTools.h" // gRefDBase
#include "pzcompel.h"           // TPZCompEl
#include "TPZRefPatternDataBase.h"

/* Command line arguments */
clarg::argInt  nsub ("-nsub", "number of substructures", 15);
clarg::argInt  porder ("-p", "porder", 2);
clarg::argInt  gcase ("-g", "gcase", 1);
clarg::argBool help ("-h", "help message", false);

void usage (char *prg)
{
    std::cout << "\nUsage: " << prg << std::endl;
    std::cout << "Arguments: "<< std::endl;
    clarg::arguments_descriptions (std::cout, "   ", "\n");
}

#ifdef PZ_LOG
#include "pzlog.h"
static TPZLogger logger("pz.Cedric-Perf");
#endif

void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial=false, const int matidtodivided=1);

RunStatsTable cedric_rst ("-ced_rdt", "TCedric::Run statistics raw data table.");


int main (int argc, char **argv)
{
    /* Parse the arguments */
    if (clarg::parse_arguments(argc, argv)) {
        std::cerr << "Error when parsing the arguments!" << std::endl;
        return 1;
    }
    /* Help message */
    if (help.get_value()) {
        usage(argv[0]);
        return 0;
    }
     /* Check if the command line arguments are valid */
    if(porder.get_value() > 10 || gcase.get_value() > 4 || nsub.get_value() > 100) {
        std::cout << "\nParameter out of avaliable limit.";
        return 1;
    }
    /* Initializing a ref patterns */
    gRefDBase.InitializeAllUniformRefPatterns();
    /* Setting interpolation order */
    TPZCompEl::SetgOrder(porder.get_value());
    /* Output file to store errors */
    std::ofstream arq("Errors.txt", std::ios::app);
    /* Create Cedric test instance */
    TCedricTest cedric;
    /* CedridTest stats */
    cedric_rst.start();
    /* Set parameters and run the test */
    cedric.Run(nsub.get_value(), gcase.get_value(), porder.get_value(), 1, arq);
    /* CedridTest stats */
    cedric_rst.stop();
}
void UniformRefinement(const int nDiv, TPZGeoMesh *gmesh, const int dim, bool allmaterial, const int matidtodivided) {
    TPZManVector<TPZGeoEl*> filhos;
    for(int D=0; D<nDiv; D++)
    {
        int nels = gmesh->NElements();
        for(int elem = 0; elem < nels; elem++)
        {
            TPZGeoEl * gel = gmesh->ElementVec()[elem];
            if(!gel || gel->HasSubElement())
                continue;
            if(dim > 0 && gel->Dimension() != dim) continue;
            if(!allmaterial) {
                if(gel->MaterialId() == matidtodivided){
                    gel->Divide(filhos);
                }
            }
            else{
                gel->Divide(filhos);
            }
#ifdef PZDEBUG
                REAL volgel = fabs(gel->Volume());
                REAL sumvol = 0.;
                for(int nsubs=0;nsubs<gel->NSubElements();nsubs++)
                    sumvol += fabs(filhos[nsubs]->Volume());
                if(!IsZero(volgel-sumvol)) {
                    std::cout << "Division of geometric element " << elem << " is wrong.\n";
                    DebugStop();
                }
#endif
        }
    }
    gmesh->ResetConnectivities();
    gmesh->BuildConnectivity();
}
