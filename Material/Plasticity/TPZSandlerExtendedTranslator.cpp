/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TPZSandlerExtendedTranslator.cpp
 * Author: thiago
 * 
 * Created on 12 de Março de 2018, 15:39
 */

#include "TPZSandlerExtendedTranslator.h"
#include "TPZChunkInTranslation.h"

TPZSandlerExtendedTranslator::TPZSandlerExtendedTranslator() {
}

TPZSandlerExtendedTranslator::TPZSandlerExtendedTranslator(const TPZSandlerExtendedTranslator& orig) {
}

void TPZSandlerExtendedTranslator::UpdateStream(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    auto old_version = chunk.mOldVersion["NeoPZ"];
    auto new_version = toVersion.at("NeoPZ");
    switch (old_version) {
        case 1:
            if (new_version <= 2){
                UpdateAttributesV1(chunk, toVersion);
                break;
            } else {
                DebugStop();
            }
        case 2:
            if (new_version == 2){
                UpdateAttributesV2(chunk, toVersion);
            } else {
                UpdateFromV2(chunk, toVersion);
            }
            break;
        default:
            UpdateAttributes(chunk, toVersion);
            break;
    }
}

void TPZSandlerExtendedTranslator::UpdateAttributesV1(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion) {
    STATE ftol;
    chunk.mOldStream.Read(&ftol);
    chunk.mNewStream.Write(&ftol);
    STATE fA;
    chunk.mOldStream.Read(&fA);
    chunk.mNewStream.Write(&fA);
    STATE fB;
    chunk.mOldStream.Read(&fB);
    chunk.mNewStream.Write(&fB);
    STATE fC;
    chunk.mOldStream.Read(&fC);
    chunk.mNewStream.Write(&fC);
    STATE fD;
    chunk.mOldStream.Read(&fD);
    chunk.mNewStream.Write(&fD);
    STATE fW;
    chunk.mOldStream.Read(&fW);
    chunk.mNewStream.Write(&fW);
    STATE fK;
    chunk.mOldStream.Read(&fK);
    chunk.mNewStream.Write(&fK);
    STATE fR;
    chunk.mOldStream.Read(&fR);
    chunk.mNewStream.Write(&fR);
    STATE fG;
    chunk.mOldStream.Read(&fG);
    chunk.mNewStream.Write(&fG);
    STATE fPhi;
    chunk.mOldStream.Read(&fPhi);
    chunk.mNewStream.Write(&fPhi);
    STATE fN;
    chunk.mOldStream.Read(&fN);
    chunk.mNewStream.Write(&fN);
    STATE fPsi;
    chunk.mOldStream.Read(&fPsi);
    chunk.mNewStream.Write(&fPsi);
    STATE fE;
    chunk.mOldStream.Read(&fE);
    chunk.mNewStream.Write(&fE);
    STATE fnu;
    chunk.mOldStream.Read(&fnu);
    chunk.mNewStream.Write(&fnu);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

void TPZSandlerExtendedTranslator::UpdateFromV2(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    STATE ftol;
    chunk.mOldStream.Read(&ftol);
    chunk.mNewStream.Write(&ftol);
    STATE fA;
    chunk.mOldStream.Read(&fA);
    chunk.mNewStream.Write(&fA);
    STATE fB;
    chunk.mOldStream.Read(&fB);
    chunk.mNewStream.Write(&fB);
    STATE fC;
    chunk.mOldStream.Read(&fC);
    chunk.mNewStream.Write(&fC);
    STATE fD;
    chunk.mOldStream.Read(&fD);
    chunk.mNewStream.Write(&fD);
    STATE fW;
    chunk.mOldStream.Read(&fW);
    chunk.mNewStream.Write(&fW);
    STATE fK;
    chunk.mOldStream.Read(&fK);
    chunk.mNewStream.Write(&fK);
    STATE fR;
    chunk.mOldStream.Read(&fR);
    chunk.mNewStream.Write(&fR);
    STATE fG;
    chunk.mOldStream.Read(&fG);
    chunk.mNewStream.Write(&fG);
    STATE fPhi;
    chunk.mOldStream.Read(&fPhi);
    chunk.mNewStream.Write(&fPhi);
    STATE fN;
    chunk.mOldStream.Read(&fN);
    chunk.mNewStream.Write(&fN);
    STATE fPsi;
    chunk.mOldStream.Read(&fPsi);
    chunk.mNewStream.Write(&fPsi);
    STATE fE;
    chunk.mOldStream.Read(&fE);
    chunk.mNewStream.Write(&fE);
    STATE fnu;
    chunk.mOldStream.Read(&fnu);
    chunk.mNewStream.Write(&fnu);
    STATE fkappa_0 = 0;
    chunk.mNewStream.Write(&fkappa_0);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

void TPZSandlerExtendedTranslator::UpdateAttributesV2(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    UpdateAttributesV1(chunk, toVersion);
}

void TPZSandlerExtendedTranslator::UpdateAttributes(TPZChunkInTranslation& chunk, const std::map<std::string, uint64_t>& toVersion){
    
    STATE ftol;
    chunk.mOldStream.Read(&ftol);
    chunk.mNewStream.Write(&ftol);
    STATE fA;
    chunk.mOldStream.Read(&fA);
    chunk.mNewStream.Write(&fA);
    STATE fB;
    chunk.mOldStream.Read(&fB);
    chunk.mNewStream.Write(&fB);
    STATE fC;
    chunk.mOldStream.Read(&fC);
    chunk.mNewStream.Write(&fC);
    STATE fD;
    chunk.mOldStream.Read(&fD);
    chunk.mNewStream.Write(&fD);
    STATE fW;
    chunk.mOldStream.Read(&fW);
    chunk.mNewStream.Write(&fW);
    STATE fK;
    chunk.mOldStream.Read(&fK);
    chunk.mNewStream.Write(&fK);
    STATE fR;
    chunk.mOldStream.Read(&fR);
    chunk.mNewStream.Write(&fR);
    STATE fG;
    chunk.mOldStream.Read(&fG);
    chunk.mNewStream.Write(&fG);
    STATE fPhi;
    chunk.mOldStream.Read(&fPhi);
    chunk.mNewStream.Write(&fPhi);
    STATE fN;
    chunk.mOldStream.Read(&fN);
    chunk.mNewStream.Write(&fN);
    STATE fPsi;
    chunk.mOldStream.Read(&fPsi);
    chunk.mNewStream.Write(&fPsi);
    STATE fE;
    chunk.mOldStream.Read(&fE);
    chunk.mNewStream.Write(&fE);
    STATE fnu;
    chunk.mOldStream.Read(&fnu);
    chunk.mNewStream.Write(&fnu);
    STATE fkappa_0;
    chunk.mOldStream.Read(&fkappa_0);
    chunk.mNewStream.Write(&fkappa_0);
    tpzElasticResponseTranslator.UpdateStream(chunk, toVersion);
}

TPZSandlerExtendedTranslator::~TPZSandlerExtendedTranslator() {
}
