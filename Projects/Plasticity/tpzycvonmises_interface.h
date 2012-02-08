// $Id: tpzycvonmises_interface.h,v 1.2 2008-03-08 03:12:53 erick Exp $
#ifndef __TPZYCVONMISES_INTERFACE_H
#define __TPZYCVONMISES_INTERFACE_H

class TPZYCVonMisesInterface
{
public:
    TPZYCVonMisesInterface() {}
    virtual ~TPZYCVonMisesInterface() {}


private:
    TPZYCVonMisesInterface( const TPZYCVonMisesInterface& source );
    void operator = ( const TPZYCVonMisesInterface& source );
};


#endif // __TPZYCVONMISES_INTERFACE_H
