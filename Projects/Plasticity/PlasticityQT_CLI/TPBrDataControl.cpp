#include "TPBrDataControl.h"

TPBrDataControl::TPBrDataControl()
{
}

int TPBrDataControl::OpenLabFile(const std::string &filename){
    std::ifstream input(filename.c_str());
    if (!input) {
        DebugStop();
    }
    int numlines = 0;
    char buf[1024];
    input.getline(buf , 1024);
    STATE x, sig_ax_t, tempo, sig_ax_dev, sig_r, eps_ax, eps_r, eps_v;
    TPBrLaboratoryData Med_txt;
    while (input) {
        input >> x >> sig_ax_t >> tempo >> sig_ax_dev >> tempo >> sig_r >> tempo >> eps_ax >> tempo >> eps_r >> tempo >> eps_v;
        if (!input) {
            break;
        }
        if(numlines >= Med_txt.fEps_Ax.size())
        {
            Med_txt.fEps_Ax.Resize(numlines+100);
            Med_txt.fSig_Ax.Resize(numlines+100);
            Med_txt.fEps_Lat.Resize(numlines+100);
            Med_txt.fSig_Lat.Resize(numlines+100);
        }
        Med_txt.fEps_Ax[numlines] = -eps_ax/100.;
        Med_txt.fSig_Ax[numlines] = -sig_ax_t;
        Med_txt.fEps_Lat[numlines] = -eps_r/100.;
        Med_txt.fSig_Lat[numlines] = -sig_r;

//    cout << "i= " << numlines << " " <<  fStressRZInput(numlines,1) << " " <<  fStressRZInput(numlines,0)<< " " << fStrainRZInput(numlines,1)<< " "<< fStrainRZInput(numlines,0) << endl;

        numlines++;
    }
    Med_txt.fEps_Ax.Resize(numlines);
    Med_txt.fSig_Ax.Resize(numlines);
    Med_txt.fEps_Lat.Resize(numlines);
    Med_txt.fSig_Lat.Resize(numlines);

    Med_txt.Set_start_idx(0);
    Med_txt.Set_end_idx(numlines);

    int pos = Medicoes.size()+1;
    Medicoes.Resize(pos);
    
    pos = pos - 1;
    Medicoes[pos] = Med_txt;
    
    //return inserted position
    return pos;
}

TPBrDataControl DADOS;
