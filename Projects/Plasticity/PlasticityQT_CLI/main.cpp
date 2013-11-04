#include "TPBrDataControl.h"
#include <iostream>

int main(int argc, char *argv[])
{
    DADOS.OpenLabFile("/home/felps/ensaio_all_columns.txt");
    std::cout << "Medicoes.NElements: " << DADOS.Medicoes.NElements();

    std::cout << DADOS.Medicoes[0].Get_start_idx();

    std::cout << "END";

}
