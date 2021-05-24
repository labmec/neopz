#include "TPZSimpleTimer.h"
#include <iostream>
#include <sstream>
void TPZSimpleTimer::Init()
{
    if(gLastTimerCreated){
        fParentTimer = gLastTimerCreated;
    }
    gLastTimerCreated = this;
    fBegin = std::chrono::high_resolution_clock::now();
}

TPZSimpleTimer::TPZSimpleTimer()
{
    Init();
}
TPZSimpleTimer::TPZSimpleTimer(const std::string name, bool alwaysPrint) :
    fName(name), fAlwaysPrint(alwaysPrint)
{
    Init();
}

TPZSimpleTimer::~TPZSimpleTimer()
{
    fEnd = std::chrono::high_resolution_clock::now();
    if(gLastTimerCreated == this){
        gLastTimerCreated = fParentTimer;
    }
    std::chrono::duration<double, std::milli> duration = fEnd-fBegin;
    std::ostringstream info{"["+fName+":"+
        std::to_string(duration.count()) +"ms"+
    fNestedTimersInfo + "]"};
    
    if(fParentTimer){
        fParentTimer->fNestedTimersInfo.append("\n\t");
        auto parent = fParentTimer;
        while(parent->fParentTimer){
            fParentTimer->fNestedTimersInfo.append("\t");
            parent = parent->fParentTimer;
        }
        fParentTimer->fNestedTimersInfo.append(info.str());
        if(fAlwaysPrint) std::cout<<info.str()<<std::endl;
    }else{
        std::cout<<info.str()<<std::endl;
    }
    
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