#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzshapepiram.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzgeoel.h"
#include "pzgmesh.h"
#include "pzstack.h"
#include "pzeltype.h"

using namespace pzshape;
using namespace std;

namespace pzrefine {


static int nsubeldata[19] = {1,1,1,1,1,3,3,3,3,3,3,3,3,9,7,7,7,7,27};

static int subeldata[19][27][2] = {//CASO DIFERENTE TAMANHO
/*00*/{{0,0}},
/*01*/{{1,1}},
/*02*/{{2,2}},
/*03*/{{3,3}},
/*04*/{{4,4}},
/*05*/{{0,5},{0,1},{1,5}},
/*06*/{{1,6},{1,2},{2,6}},
/*07*/{{2,7},{2,3},{3,7}},
/*08*/{{0,8},{0,3},{3,8}},
/*09*/{{0,9},{0,4},{4,9}},
/*10*/{{1,10},{1,4},{4,10}},
/*11*/{{2,11},{2,4},{4,11}},
/*12*/{{3,12},{3,4},{4,12}},
// /*13*/{{0,13},{1,13},{2,13},{3,13},{0,2},{0,6},{0,7},{1,7},{2,8}},
// /*14*/{{0,14},{1,14},{4,14},{6,11},{0,10},{1,9},{4,5}},
// /*15*/{{1,15},{2,15},{4,15},{7,10},{1,11},{2,10},{4,6}},
// /*16*/{{2,16},{3,16},{4,16},{8,13},{2,12},{3,11},{4,7}},
// /*17*/{{0,17},{3,17},{4,17},{9,12},{0,12},{3,9},{4,8}},
/*13*/{{0,2},{0,6},{0,7},{1,7},{2,8},{0,13},{1,13},{2,13},{3,13}},
/*14*/{{0,10},{1,9},{4,5},{0,14},{1,14},{4,14},{6,11}},
/*15*/{{1,11},{2,10},{4,6},{7,10},{1,15},{2,15},{4,15}},
/*16*/{{2,12},{3,11},{4,7},{2,16},{3,16},{4,16},{8,13}},
/*17*/{{0,12},{3,9},{4,8},{0,17},{3,17},{4,17},{9,12}},
/*18*/{{0,18},{1,18},{2,18},{3,18},{4,18},{5,18},{6,14},{7,14},{8,14},
       {9,14},{0,11},{1,12},{2,9},{3,10},{0,15},{0,16},{1,16},{1,17},
       {2,14},{2,17},{3,14},{3,15},{4,13},{6,12},{7,13},{8,11},{9,10}}
};


static int MidSideNodes[9][2]  = { 
	{0,1},{1,2},{2,3},
	{3,0},{4,0},{4,1},
	{4,2},{4,3},{0,2} 
};

static REAL MidCoord[9][3] = { 
	{0.,-1.,0.},{1.,0.,0.},{0.,1.,0.},
	{-1.,0.,0.},{-.5,-.5,.5},{.5,-.5,.5},
	{.5,.5,.5},{-.5,.5,.5},{0.,0.,0.} 
};

/**
 * define as conectividades entre sub-elementos
 * linha i � filho i, {a,b,c} = {lado do filho atual,
 * irm�o vizinho,lado do vizinho}
 */
const int NumInNeigh = 18;
static int InNeigh[10][NumInNeigh][3] = { 
	{{1,6,1},{2,6,2},{3,9,3},{4,6,0},{6,6,5},{7,9,7},{10,6,4},{11,6,6},{12,9,8},{15,6,10},{16,9,11},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
	{{0,0,1},{2,7,0},{3,7,3},{4,7,1},{7,7,7},{8,0,6},{9,6,8},{11,7,4},{12,7,8},{16,7,11},{17,6,12},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
	{{0,8,1},{1,1,2},{3,8,2},{4,8,3},{5,1,7},{8,8,5},{9,8,8},{10,7,6},{12,8,9},{14,7,13},{17,8,12},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,0,3},{1,9,0},{2,2,3},{4,9,2},{5,0,7},{6,2,8},{9,9,9},{10,9,6},{11,8,6},{14,9,13},{15,8,10},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
	{{0,0,4},{1,6,3},{2,5,3},{3,3,4},{5,6,7},{6,7,5},{7,8,7},{8,9,5},{13,5,13},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,4,1},{1,9,1},{2,8,0},{3,7,2},{4,0,2},{5,4,5},{6,4,8},{7,4,7},{8,4,6},{9,6,9},{10,9,4},{11,8,4},{12,7,9},{13,4,13},{14,6,13},{15,9,10},{16,8,11},{17,7,12}},
    {{0,5,1},{1,1,0},{2,1,3},{3,1,4},{4,0,10},{5,1,8},{6,5,10},{7,5,5},{8,1,9},{9,1,12},{10,0,15},{12,1,17},{13,5,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,2,1},{1,5,0},{2,2,4},{3,2,0},{4,1,11},{5,5,8},{6,2,10},{7,2,5},{8,5,9},{9,2,9},{11,1,16},{12,5,17},{13,2,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,4,3},{1,3,1},{2,3,2},{3,4,2},{4,3,10},{5,3,6},{6,3,11},{7,5,7},{8,5,12},{9,2,12},{10,3,15},{11,5,16},{12,2,17},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}},
    {{0,5,4},{1,4,0},{2,5,2},{3,3,0},{4,0,11},{5,5,6},{6,5,11},{7,3,5},{8,0,12},{9,3,9},{10,5,15},{11,0,16},{13,3,14},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1},{-1,-1,-1}} 
};



/**
 * define os cantos locais dos fihos
 */
static int CornerSons[10][5] = { 
	{0,5,13,8,9},{5,1,6,13,10},
	{13,6,2,7,11},{8,13,7,3,12},
	{9,10,11,12,4},{10,9,12,11,13},
	{9,5,13,10,-1},{6,10,11,13,-1},
	{12,13,7,11,-1},{13,9,12,8,-1} 
};

static REAL buildt[10][19][4][3] = {//por colunas
/*S0*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,-1,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*06*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*10*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*11*/{{-0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*12*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,-.5,0.}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.,1,0.},{-0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*16*/{{1,0.,0.},{0.5,-0.5,0.5},{0.,0.,0.},{-1,0.,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,-.5,0.}}},
/*S1*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,-1,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*11*/{{-.25,.25,0.},{0.,0.,0.},{0.,0.,0.},{.25,.25,0.}},
      /*12*/{{0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,-.5,0.}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{1,0.,0.},{0.5,-0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,-.5,0.}}},
/*S2*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{1,1,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*06*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.,-.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*09*/{{0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*10*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*12*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{.5,.5,0.}},
      /*14*/{{1,0.,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*17*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{.5,.5,0.}}},
/*S3*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{-1,1,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*06*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*07*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{-0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*11*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*13*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,0.},{-.5,.5,0.}},
      /*14*/{{1,0.,0.},{0.5,0.5,0.5},{0.,0.,0.},{-1,0.,0.}},
      /*15*/{{0.,1,0.},{-0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{-.5,.5,0.}}},
/*S4*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,1}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*09*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*10*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*11*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*12*/{{.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*13*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*14*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*15*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*16*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*17*/{{0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*18*/{{.5,0.,0.},{0.,.5,0.},{0.,0.,.5},{0.,0.,.5}}},
/*S5*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*05*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{-0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*09*/{{-0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*10*/{{0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*11*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*12*/{{-0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*13*/{{-0.5,0.,0.},{0.,0.5,0.},{0.,0.,0.},{0.,0.,0.5}},
      /*14*/{{-1,0.,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*15*/{{0.,1,0.},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*16*/{{-1,0.,0.},{-0.5,-0.5,-0.5},{0.,0.,0.},{0.5,0.5,0.5}},
      /*17*/{{0.,1,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*18*/{{-.5,0.,0.},{0.,.5,0.},{0.,0.,-.5},{0.,0.,.5}}},
/*S6*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*05*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,-.5,0.}},
      /*06*/{{-0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*09*/{{0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*10*/{{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*11*/{{0.5,-0.5,0.},{0.5,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*12*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,-1,0.}},
      /*13*/{{0.5,0.5,-0.5},{1,0.,0.},{0.,0.,0.},{-0.5,-0.5,0.5}},
      /*14*/{{.5,-.5,-.5},{.5,.5,-.5},{1.,.0,.0},{-.5,-.5,.5}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S7*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{.5,0.,0.}},
      /*08*/{{-0.25,0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,-0.25,0.25}},
      /*09*/{{-0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*10*/{{-0.5,0.5,0.},{0.,0.5,0.},{0.,0.,0.},{0.5,0.,0.}},
      /*11*/{{-0.5,-0.5,0.5},{-1,0.,0.},{0.,0.,0.},{1,0.,0.}},
      /*12*/{{0.,1,0.},{-0.5,0.5,-0.5},{0.,0.,0.},{0.5,-0.5,0.5}},
      /*13*/{{-0.5,0.5,0.5},{-1,0.,0.},{0.,0.,0.},{1,0.,0.}},
      /*14*/{{-.5,-.5,.5},{-.5,.5,.5},{-1.,0.,0.},{1.,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S8*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*05*/{{0.,.5,0.},{0.,0.,0.},{0.,0.,0.},{0.,.5,0.}},
      /*06*/{{-0.25,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*07*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*08*/{{0.25,0.25,0.25},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.25}},
      /*09*/{{0.,0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.,0.,0.},{-0.5,0.5,0.5}},
      /*11*/{{0.5,-0.5,-0.5},{1,0.,0.},{0.,0.,0.},{-0.5,0.5,0.5}},
      /*12*/{{0.,1,0.},{0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*13*/{{0.5,-0.5,0.},{0.5,0.,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*14*/{{.5,-.5,-.5},{.5,.5,-.5},{1.,0.,0.},{-.5,.5,.5}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}},
/*S9*/{
      /*00*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*01*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*02*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*03*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*04*/{{-0.25,-0.25,0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,-0.25,0.25}},
      /*05*/{{0.25,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.5,0.}},
      /*06*/{{0.25,-0.25,-0.25},{0.,0.,0.},{0.,0.,0.},{-0.25,0.25,0.25}},
      /*07*/{{-.5,0.,0.},{0.,0.,0.},{0.,0.,0.},{-.5,0.,0.}},
      /*08*/{{0.25,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.25,0.25,0.}},
      /*09*/{{0.,-0.25,0.},{0.,0.,0.},{0.,0.,0.},{0.5,0.25,0.}},
      /*10*/{{-0.5,-0.5,0.5},{-0.5,0.5,0.5},{0.,0.,0.},{0.,0.,0.}},
      /*11*/{{-0.5,-0.5,0.5},{-1,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*12*/{{0.5,0.,0.},{0.5,-0.5,0.},{0.,0.,0.},{0.,0.5,0.}},
      /*13*/{{-0.5,0.5,0.5},{-1,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*14*/{{-.5,-.5,.5},{-.5,.5,.5},{-1.,0.,0.},{0.,0.,0.}},
      /*15*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*16*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*17*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}},
      /*18*/{{0.,0.,0.},{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}}}
};

//Into Divides is necesary to consider the connectivity with the all neighboards
void TPZRefPyramid::Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec) {
	int i;
	if(geo->HasSubElement()) {
		SubElVec.Resize(NSubEl);
		for(i=0;i<NSubEl;i++) SubElVec[i] = geo->SubElement(i);
		return;//If exist fSubEl return this sons
	}
	int j,sub,matid=geo->MaterialId(),index;
	int np[TPZShapePiram::NSides];//guarda conectividades dos 8 subelementos
	for(j=0;j<TPZShapePiram::NCornerNodes;j++) np[j] = geo->NodeIndex(j);
	for(sub=TPZShapePiram::NCornerNodes;sub<14;sub++) {
		NewMidSideNode(geo,sub,index);
		np[sub] = index;
	}
	// creating new subelements
	for (i=0;i<6;i++){
	  TPZManVector<int> cornerindexes(TPZShapePiram::NCornerNodes);
	  for(int j=0;j<TPZShapePiram::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
	  TPZGeoEl *pi3sub = geo->Mesh()->CreateGeoElement(EPiramide,cornerindexes,matid,index);
	  geo->SetSubElement(i,pi3sub);
	}
	for (;i<10;i++){
	  TPZManVector<int> cornerindexes(TPZShapeTetra::NCornerNodes);
	  for(int j=0;j<TPZShapeTetra::NCornerNodes;j++) cornerindexes[j] = np[CornerSons[i][j]];
	  TPZGeoEl *t3sub = geo->Mesh()->CreateGeoElement(ETetraedro,cornerindexes,matid,index);
	  geo->SetSubElement(i,t3sub);
	}
	SubElVec.Resize(NSubEl);
	for(sub=0;sub<NSubEl;sub++) {
		SubElVec[sub] = geo->SubElement(sub);
		SubElVec[sub]->SetFather(geo);
		SubElVec[sub]->SetFather(geo->Index());
	}
	for(i=0;i<NSubEl;i++) {//conectividades entre os filhos : viz interna
		for(j=0;j<NumInNeigh;j++) {        //lado do subel  numero do filho viz.             lado do viz.
      		int elside = InNeigh[i][j][0];//lado do subel
			if(elside == -1) break;
 			geo->SubElement(i)->SetNeighbour(elside,TPZGeoElSide(geo->SubElement(InNeigh[i][j][1]),InNeigh[i][j][2]));
		}
	}
	geo->SetSubElementConnectivities();
}

void TPZRefPyramid::NewMidSideNode(TPZGeoEl *gel,int side,int &index) {

	MidSideNodeIndex(gel,side,index);
	if(index < 0) {
		TPZGeoElSide gelside = gel->Neighbour(side);
		if(gelside.Element()) {
			while(gelside.Element() != gel) {
				gelside.Element()->MidSideNodeIndex(gelside.Side(),index);
				if(index!=-1) return;
				gelside = gelside.Neighbour();
			}	
		}
		TPZVec<REAL> par(3,0.);
		TPZVec<REAL> coord(3,0.);
		if(side < TPZShapePiram::NCornerNodes) {
			index = gel->NodeIndex(side); 
			return;
		}
		//aqui side = 8 a 26
		side-=TPZShapePiram::NCornerNodes;//0,1,..,18
		par[0] = MidCoord[side][0];
		par[1] = MidCoord[side][1];
		par[2] = MidCoord[side][2];
		gel->X(par,coord);
		index = gel->Mesh()->NodeVec().AllocateNewElement();
		gel->Mesh()->NodeVec()[index].Initialize(coord,*gel->Mesh());
	}
}

void TPZRefPyramid::MidSideNodeIndex(TPZGeoEl *gel,int side,int &index) {
	index = -1;
	if(side<0 || side>TPZShapePiram::NSides-1) {
		PZError << "TPZRefPyramid::MidSideNodeIndex. Bad parameter side = " << side << endl;
		return;
	}
	//sides 0 a 7
	if(side<TPZShapePiram::NCornerNodes) {//o n� medio do lado 0 � o 0 etc.
		index = (gel)->NodeIndex(side);
		return; 
	}
	//o n� medio da face � o centro e o n� medio do centro � o centro
	//como n� de algum filho se este existir
	//caso tenha filhos � o canto de algum filho, se n�o tiver filhos retorna -1
	if(gel->HasSubElement()) {
		side-=TPZShapePiram::NCornerNodes;
		if(side >= NSubEl) {
			index = -1;
			PZError << "TPZRefPyramid : MidSideNodeIndex called for wrong side\n";
		}
		index=(gel->SubElement(MidSideNodes[side][0]))->NodeIndex(MidSideNodes[side][1]);
	}
}

void TPZRefPyramid::GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel){

	subel.Resize(0);
	if(side<0 || side>TPZShapePiram::NSides || !father->HasSubElement()){
		PZError << "TPZRefPyramid::GetSubelements2 called with error arguments\n";
		return;
	}
	int nsub = NSideSubElements(side);//nsubeldata[side];
	for(int i=0;i<nsub;i++)
		subel.Push(TPZGeoElSide(father->SubElement(subeldata[side][i][0]),
												subeldata[side][i][1]));
}

int TPZRefPyramid::NSideSubElements(int side) {  
	if(side<0 || side>TPZShapePiram::NSides-1){
		PZError << "TPZRefPyramid::NSideSubelements2 called with error arguments\n";
		return -1;
	}
	return nsubeldata[side];
}

//int TPZRefPyramid::NSideSubElements(int side) {
//  if(side < 0 || side > 26) {
//    PZError << "TPZRefPyramid::NSideSubElements called for side " << side << endl;
//    return 0;
//  }
//  if(side==26) return 8;//centro
//  if(side>19 && side<26) return 4;//faces
//  if(side>7) return 2;//lados
//  return 1;//cantos
//}


TPZTransform TPZRefPyramid::GetTransform(int side,int whichsubel){
	if(side<0 || side>TPZShapePiram::NSides-1){
		PZError << "TPZRefPyramid::GetTransform side out of range or father null\n";
		return TPZTransform(0,0);
	}
	int smalldim;
	if(whichsubel <6) smalldim = TPZShapePiram::SideDimension(side);
	else smalldim = TPZShapeTetra::SideDimension(side);
	int fatherside = FatherSide(side,whichsubel);
	int largedim = TPZShapePiram::SideDimension(fatherside);
	TPZTransform trans(largedim,smalldim);
	int i,j;
	for(i=0; i<largedim; i++) {
		for(j=0; j<smalldim; j++) {
			trans.Mult()(i,j) = buildt[whichsubel][side][j][i];
		}
		trans.Sum() (i,0) = buildt[whichsubel][side][3][i];
	}
	return trans;
}

static int fatherside[10][19] = {
/*00*/{0,5,13,8,9,5,13,13,8,9,14,18,17,13,14,18,18,17,18},
/*01*/{5,1,6,13,10,5,6,13,13,14,10,15,18,13,14,15,18,18,18},
/*02*/{13,6,2,7,11,13,6,7,13,18,15,11,16,13,18,15,16,18,18},
/*03*/{8,13,7,3,12,13,13,7,8,17,18,16,12,13,18,18,16,17,18},
/*04*/{9,10,11,12,4,14,15,16,17,9,10,11,12,18,14,15,16,17,18},
/*05*/{10,9,12,11,13,14,17,16,15,18,18,18,18,18,18,18,18,18,18},
/*06*/{9,5,13,10,14,13,18,14,14,18,18,14,18,18,18,-1,-1,-1,-1},
/*07*/{6,10,11,13,15,15,15,13,18,18,15,18,18,18,18,-1,-1,-1,-1},
/*08*/{12,13,7,11,18,13,16,16,18,16,18,18,18,16,18,-1,-1,-1,-1},
/*09*/{13,9,12,8,18,17,18,13,17,17,18,18,17,18,18,-1,-1,-1,-1},
};

// static int fatherside2[2][19] = {//pir�mide filho de tetraedro
  // /*04*/{4,8,9,6,7,11,12,13,10,11,11,13,13,14,11,14,13,14,14},
  // /*05*/{8,4,6,9,5,11,10,13,12,12,10,10,12,14,14,10,14,12,14} };

int TPZRefPyramid::FatherSide(int side,int whichsubel){


	if(side<0 || side>TPZShapePiram::NSides-1 || whichsubel <0 || whichsubel >= TPZRefPyramid::NSubEl){
		PZError << "TPZRefPyramid::Father2 called error" << endl;
		return -1;
	}
	return fatherside[whichsubel][side];
}

};
