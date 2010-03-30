//$Id: TPZGuiInterface.cpp,v 1.1 2010-03-30 14:37:06 fortiago Exp $

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



