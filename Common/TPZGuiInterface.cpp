//$Id: TPZGuiInterface.cpp,v 1.2 2010-04-06 17:22:03 fortiago Exp $

#include "TPZGuiInterface.h"
#include "pzerror.h"

TPZGuiInterface::TPZGuiInterface(){
	this->fCanceled = false;
	this->fMessage = "";
	this->fProgressBarPos = 0;
	this->fProgressBarMaxPos = 0;
	this->fProgressBarMinPos = 0;
}

TPZGuiInterface::~TPZGuiInterface(){
	///nothing to be done
	int i = 0;
}

void TPZGuiInterface::UpdateCaption(){
	std::cout << fMessage.c_str() << "\n"
			  << "Progress bar = " << fProgressBarPos << "/" << fProgressBarMaxPos
			  << "\n";
}

void TPZGuiInterface::Start(){
	std::cout << "Starting execution\n";
}

void TPZGuiInterface::End(){
	std::cout << "Execution finished\n";
}

void TPZGuiInterface::ShowErrorMessage(std::string message){
	PZError << message.c_str() << "\n";
	DebugStop();
}


void TPZGuiInterface::SetKilled(){
	this->fCanceled = true;
}


bool TPZGuiInterface::AmIKilled(){
	return this->fCanceled;
}



