#include "TPZSimpleTimer.h"
#include <iostream>

TPZSimpleTimer::TPZSimpleTimer() : fBegin(std::chrono::high_resolution_clock::now())
{
}

TPZSimpleTimer::TPZSimpleTimer(const std::string name) : fName(name), fBegin(std::chrono::high_resolution_clock::now())
{
}

TPZSimpleTimer::~TPZSimpleTimer()
{
    fEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = fEnd-fBegin;
    std::cout<<"["<<fName<<","<<duration.count()<<"]"<<std::endl;
}

double TPZSimpleTimer::ReturnTimeDouble()
{
    auto now = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = now-fBegin;
    return duration.count();
}

std::string TPZSimpleTimer::ReturnTimeString()
{
    return std::to_string(ReturnTimeDouble());
}