//
// C++ Interface: pzmaterialid
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@corona>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZMATERIALIDH
#define PZMATERIALIDH

const int TPZMATERIALID = 300;
const int TPZDISCONTINUOUSGALERKIN = 301;
const int TPZMAT2DLINID = 302;
const int TPZCONSERVATIONLAW2ID = 303;
const int TPZEULERCONSLAW2ID = 304;
const int TPZARTDIFFID = 305;
const int TPZBNDCONDID = 306;
const int TPZELASTICITYMATERIALID = 307;
const int TPZELASTICITY3DMATERIALID = 308;
const int TPZMATTEST3DID = 309;
const int TPZNLMAT1D = 310;
const int TPZMATPOISSON3D = 311;



void RegisterMaterialClasses();


#endif //PZMATERIALIDH
