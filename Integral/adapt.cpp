
#include "adapt.h"


const REAL Adapt::alpha = sqrt(2.0/3.0);
const REAL Adapt::beta = 1.0/sqrt(5.0);
const REAL Adapt::x1 = 0.942882415695480;
const REAL Adapt::x2 = 0.641853342345781;
const REAL Adapt::x3 = 0.236383199662150;
const REAL Adapt::x[12] = {0, -x1, -alpha, -x2, -beta, -x3, 0.0, x3, beta, x2, alpha, x1};
