/* 
 * File:   TPZChunkTranslator.h
 * Author: thiago
 *
 * Created on 7 de Março de 2018, 17:14
 */

#ifndef TPZCHUNKTRANSLATOR_H
#define TPZCHUNKTRANSLATOR_H

#include <map>
#include <string>

class TPZChunkInTranslation;

class TPZChunkTranslator {
  
    
public:
  
    virtual ~TPZChunkTranslator();
  
    virtual void UpdateStream(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion) {
        this->UpdateAttributes(chunk, toVersion);
    }
    
    virtual void UpdateAttributes(TPZChunkInTranslation &chunk, const std::map<std::string, uint64_t> &toVersion)=0;
    
    virtual void SetClassId(int classid){
        TPZChunkTranslator::classid = classid;
    }
    
    virtual int GetClassId() const {
        return classid;
    }
    
private :
    static int classid;
};

#endif /* TPZCHUNKTRANSLATOR_H */

