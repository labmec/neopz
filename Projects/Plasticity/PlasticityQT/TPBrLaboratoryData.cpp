#include "TPBrLaboratoryData.h"
#include "TPZPlasticitySimulation.h"
#include "TPBrDataControl.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("QT.laboratorydata"));
#endif

TPBrLaboratoryData::TPBrLaboratoryData()
{
    fstart_idx = -1;
    fend_idx = -1;
    felastic_idx = -1;
}

TPBrLaboratoryData::TPBrLaboratoryData(const std::string &filename)
{
    ReadInputStrainStress(filename);
    //setting start and end idxs to first and end points
    fstart_idx = 0;
    fend_idx = this->fSig_Ax.size() - 1;
    felastic_idx = 0;
}

int TPBrLaboratoryData::InsertSimulationData (const TPBrSimulationData &simdataobj) {
      int simid = DADOS.GenerateNewIndex();
      fSimulacoes[simid] = simdataobj;
      fSimulacoes[simid].SetGlobalId(simid);
      return simid;
}

int TPBrLaboratoryData::RunSimulation (TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj) {
  
  // RUN SIMULATION AND RETURN ITS INDEX
  
  TPZPlasticitySimulation newSimulation;
  
  newSimulation.ReadInputStrainStress(fSig_Ax, fEps_Ax, fSig_Lat, fEps_Lat);
  newSimulation.SetSimulationInitialStep(fstart_idx);
  // TRANSICAO ELASTICA AQUIIIIIIIIIIIIIIIII!!
  newSimulation.SetSandlerDimaggio(obj);
  newSimulation.SetElasticTransition(felastic_idx);
  newSimulation.PerformSimulation();
    TPZVec<REAL> sigax,sigr,epsax,epsr;
    newSimulation.GetSimulatedStrainStress(sigax, epsax, sigr, epsr);
		
		if(felastic_idx < fstart_idx)
		{
			std::cout << "The elastic index should always be larger than the pore closure index\n";
		}
		else
		{
			REAL errax = sigax[felastic_idx-fstart_idx]-fSig_Ax[felastic_idx-fstart_idx];
			REAL errlat = sigr[felastic_idx-fstart_idx]-fSig_Lat[felastic_idx-fstart_idx];
			std::cout << "errax " << errax << " errlat " << errlat << std::endl;
		}  
//#ifdef LOG4CXX
//		if(logger->isDebugEnabled())
//		{
//			std::stringstream sout;
//			sout << "epsax " << epsax << std::endl << "epsr " << epsr << std::endl << "sigax " << sigax << std::endl << "sigr " << sigr;
//			LOGPZ_DEBUG(logger,sout.str())
//		}
			
//#endif
    
    TPBrSimulationData result;
    int medid = GlobalId();
    result.Set_medicao_idx(medid);
    result.Set_start_idx(fstart_idx);
    result.Set_end_idx(fend_idx);
    result.Set_elastic_trans_idx(felastic_idx);
    result.SetStrainStress(sigax, epsax, sigr, epsr);
    int resultid = DADOS.InsertSimulationData(medid, result);
    
  //return inserted position
  return resultid;
}

/// read the input strain and stress from the laboratory file
void TPBrLaboratoryData::ReadInputStrainStress(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input) {
        DebugStop();
    }
    int numlines = 0;
    char buf[1024];
    input.getline(buf , 1024);
    STATE x, sig_ax_t, tempo, sig_ax_dev, sig_r, eps_ax, eps_r, eps_v;
    while (input) {
        input >> x >> sig_ax_t >> tempo >> sig_ax_dev >> tempo >> sig_r >> tempo >> eps_ax >> tempo >> eps_r >> tempo >> eps_v;
        if (!input) {
            break;
        }
        if(numlines >= fSig_Ax.size())
        {
            fSig_Ax.Resize(numlines+100, 2);
            fSig_Lat.Resize(numlines+100, 2);
            fEps_Ax.Resize(numlines+100, 2);
            fEps_Lat.Resize(numlines+100, 2);
        }
        fSig_Ax[numlines] = -sig_ax_t;
        fSig_Lat[numlines] = -sig_r;
        fEps_Ax[numlines] = -eps_ax/100.;
        fEps_Lat[numlines] = -eps_r/100.;
        
        //    cout << "i= " << numlines << " " <<  fStressRZInput(numlines,1) << " " <<  fStressRZInput(numlines,0)<< " " << fStrainRZInput(numlines,1)<< " "<< fStrainRZInput(numlines,0) << endl;
        
        numlines++;
    }
    fSig_Ax.Resize(numlines);
    fSig_Lat.Resize(numlines);
    fEps_Lat.resize(numlines);
    fEps_Ax.resize(numlines);
}

    void TPBrLaboratoryData::IdentifyElasticity (REAL &Young, REAL &Poisson) {
        //your code here
        //will work already in GUI

        REAL SigAx_start = this->fSig_Ax[fstart_idx];
        REAL SigLat_start = this->fSig_Lat[fstart_idx];
        REAL EpsAx_start = this->fEps_Ax[fstart_idx];
        REAL EpsLat_start = this->fEps_Lat[fstart_idx];

        REAL SigAx_elast = this->fSig_Ax[felastic_idx];
        REAL SigLat_elast = this->fSig_Lat[felastic_idx];
        REAL EpsAx_elast = this->fEps_Ax[felastic_idx];
        REAL EpsLat_elast = this->fEps_Lat[felastic_idx];
				
				REAL DSigR = SigLat_elast-SigLat_start;
				REAL DEpsR = EpsLat_elast - EpsLat_start;
				REAL DSigAx = SigAx_elast-SigAx_start;
				REAL DEpsAx = EpsAx_elast-EpsAx_start;
				
				Young = ((DSigAx - DSigR)*(DSigAx + 2*DSigR))/(-2*DEpsR*DSigR + DEpsAx*(DSigAx + DSigR));
				Poisson = (-(DEpsR*DSigAx) + DEpsAx*DSigR)/(-2*DEpsR*DSigR + DEpsAx*(DSigAx + DSigR));
				
//#ifdef DEBUG
				{
					
					TPZElasticResponse ER;
					ER.SetUp(Young,Poisson);
					TPZTensor<STATE> eps_start,eps_elast,sig_start,sig_elast,del_eps, del_sig, del_sig_input;
					eps_start.XX() = EpsLat_start;
					eps_start.YY() = EpsLat_start;
					eps_start.ZZ() = EpsAx_start;
					eps_elast.XX() = EpsLat_elast;
					eps_elast.YY() = EpsLat_elast;
					eps_elast.ZZ() = EpsAx_elast;
					sig_start.XX() = SigLat_start;
					sig_start.YY() = SigLat_start;
					sig_start.ZZ() = SigAx_start;
					sig_elast.XX() = SigLat_elast;
					sig_elast.YY() = SigLat_elast;
					sig_elast.ZZ() = SigAx_elast;
					del_eps = eps_elast;
					del_eps -= eps_start;
					del_sig_input = sig_elast;
					del_sig_input -= sig_start;
					REAL epsxx = del_sig_input.XX()/Young - del_sig_input.YY()*Poisson/Young - del_sig_input.ZZ()*Poisson/Young;
					ER.Compute(del_eps,del_sig);
					del_sig -= del_sig_input;
#ifdef LOG4CXX
					{
						std::stringstream sout;
						sout << "del_sig_input " << del_sig_input << std::endl;
						sout << "error " << del_sig << std::endl;
						LOGPZ_DEBUG(logger,sout.str())
					}
#endif
					
				}
//#endif

    }

