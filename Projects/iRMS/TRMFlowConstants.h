//---------------------------------------------------------------------------

#ifndef TRMFlowConstantsH
#define TRMFlowConstantsH


const int _ReservMatId = 1;
const int _WellMatId3D = 2;
const int _WellFacesMatId = 3;
const int _LateralReservBC = 4;
const int _ConfinementReservBC = 5;
const int _WellMatId1D = 6;
const int _PerforationMatId = 7;
const int _LinerMatId = 8;
//nos usamos _WellMatId3DIsLeft ou _WellMatId3DIsRight const int _SkinMatId = 9;///para skin e canhoneado
const int _BlendArcMatId = 10;
const int _FractureMatId = 11;
const int _ReservMatId1DFake = 12;
const int _WellMatId1DFake = 13;
const int _WellHeelMatId = 14;
const int _NullPermMatId = 1123581321;

const int _mioloReservFaces = 16;
const int _WellToeMatId = 15;
const int _WellMatId3DIsLeft = 17;
const int _WellMatId3DIsRight = 19;
const int _waterLayerMatId = 21;

//temporary materials for directional refinements tecniche in geomesh generation
const int _wellTipMat = 22;//temporary material id
const int _mioloBottomMat = 23;//temporary material id
const int _mioloTopMat = 24;//temporary material id
const int _mioloLateralMat = 25;//temporary material id
const int _mioloBoundary = 37;
const int _ConfinementReservBCbottom = 26;//temporary material id
const int _ConfinementReservBCtop = 27;//temporary material id
const int _Reservoir_regiaoMiolo = 28;//temporary material id
const int _Reservoir_extremidades1 = 29;//temporary material id
const int _Reservoir_extremidades2 = 30;//temporary material id
const int _Water_regiaoMiolo = 31;//temporary material id
const int _Water_extremidades1 = 32;//temporary material id
const int _Water_extremidades2 = 33;//temporary material id
const int _WaterMioloBottomMat = 34;//temporary material id
const int _ellipseArcMat = 35;

//Additional constants
const int _ReservoirInletPressure = 36;
const int _ReservoirOutletPressure = 37;
const int _ReservoirNonFluxBoundary = 38;

// two dimensional elements between the 3D well elements and the reservoir
const int _Well3DReservoirFaces = 39;

const int _ReservoirInterface = 40;



//---------------------------------------------------------------------------
#endif
