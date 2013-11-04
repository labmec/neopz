#include "TPBrDataControl.h"
#include <iostream>

int main(int argc, char *argv[])
{
  
    string FileName, dirname = PZSOURCEDIR;
    FileName = dirname + "/Projects/Plasticity/";
    FileName += "ensaio_all_columns.txt";
  
    DADOS.OpenLabFile(FileName);
    std::cout << "Medicoes.NElements: " << DADOS.Medicoes.NElements() << std::endl;

    std::cout << DADOS.Medicoes[0].Get_start_idx()  << std::endl;

    std::cout << "END" << std::endl;

}
