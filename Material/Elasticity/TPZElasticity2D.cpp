#include "Elasticity/TPZElasticity2D.h"
#include "pzaxestools.h"



TPZElasticity2D::TPZElasticity2D() : 
TPZRegisterClassId(&TPZElasticity2D::ClassId),
TBase(), ff(3,0.) {
	fE_def	= -1.;  // Young modulus
	fnu_def	= -1.;   // poisson coefficient
	ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
    
	
	//Added by Cesar 2001/03/16
	fPreStressXX = 0.;  //Prestress in the x direction
	fPreStressYY = 0.;  //Prestress in the y direction
	fPreStressXY = 0.;  //Prestress in the z direction
	fPreStressZZ = 0.;  //Prestress in the z direction
	fPlaneStress = 0;
    
    // Added by Philippe 2012
    fPostProcIndex = 0;
    fUseVecShape = false;
}

TPZElasticity2D::TPZElasticity2D(const TPZElasticity2D &copy) : TBase(copy), TPZMatErrorSingleSpace<STATE>(copy), fE_def(copy.fE_def), fnu_def(copy.fnu_def), fConstitutiveLaw(copy.fConstitutiveLaw),
fElasticity(copy.fElasticity), ff(copy.ff), fPreStressXX(copy.fPreStressXX),
fPreStressYY(copy.fPreStressYY), fPreStressXY(copy.fPreStressXY), fPreStressZZ(copy.fPreStressZZ),
fPlaneStress(copy.fPlaneStress), fUseVecShape(copy.fUseVecShape)
{
    if(copy.HasExactSol()) {
        SetExactSol(copy.ExactSol(), copy.PolynomialOrderExact());
    }
//    std::cout << "copy has exact " << copy.HasExactSol() << std::endl;
//    std::cout << "I have exact sol " << HasExactSol() << std::endl;
}

TPZElasticity2D::TPZElasticity2D(int id) :
    TPZRegisterClassId(&TPZElasticity2D::ClassId),
    TBase(id), ff(3,0.) {
	fE_def	= -1.;  // Young modulus
	fnu_def	= -1.;   // poisson coefficient
	ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
    
	
	//Added by Cesar 2001/03/16
	fPreStressXX = 0.;  //Prestress in the x direction
	fPreStressYY = 0.;  //Prestress in the y direction
	fPreStressXY = 0.;  //Prestress in the z direction
	fPreStressZZ = 0.;  //Prestress in the z direction
	fPlaneStress = 0;
    
    // Added by Philippe 2012
    fPostProcIndex = 0;
}


TPZElasticity2D::TPZElasticity2D(int id, STATE E, STATE nu,
                                 STATE fx, STATE fy, int planestress) : TPZRegisterClassId(&TPZElasticity2D::ClassId), TBase(id), ff(3,0.){
    // std::cout << __PRETTY_FUNCTION__ << std::endl;
    fConstitutiveLaw.SetIsotropicProperties(E, nu);
    fConstitutiveLaw.SetPlaneStress(planestress);
    fConstitutiveLaw.SetDimension(2);
    fPlaneStress = planestress;
    fE_def = E;
    fnu_def = nu;
    ff[0] = fx;
    ff[1] = fy;
    ff[2] = 0.; // Z component of the body force - not used for this class
    
    fPreStressXX = 0.;  //Prestress in the x direction
    fPreStressYY = 0.;  //Prestress in the y direction
    fPreStressXY = 0.;  //Prestress in the z direction
    fPreStressZZ = 0.;  //Prestress in the z direction
    
    // Added by Philippe 2012
    fPostProcIndex = 0;
    fUseVecShape = false;

}

void TPZElasticity2D::SetPreStress(STATE Sigxx, STATE Sigyy, STATE Sigxy, STATE Sigzz){
	fPreStressXX = Sigxx;
	fPreStressYY = Sigyy;
	fPreStressXY = Sigxy;
    fPreStressZZ = Sigzz;
}

int TPZElasticity2D::NStateVariables() const {
	return 2;
}

void TPZElasticity2D::Print(std::ostream &out) const {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
    if(fElasticity)
    {
        out << "Elasticity coeficients are determined by a function\n";
    }
    else
    {
        out << "\tE   = " << fE_def   << '\n';
        out << "\tnu   = " << fnu_def   << '\n';
    }
	out << "\tF   = " << ff[0] << ' ' << ff[1]   <<'\n';
	out << "\t PreStress: \n"
	<< "Sigma xx = \t" << fPreStressXX << "\t"
	<< "Sigma yy = \t" << fPreStressYY << "\t"
        << "Sigma xy = \t" << fPreStressXY << "Sigma zz = \t" << fPreStressZZ << std::endl;
}

/// @brief Compute the constitutive matrix D
/// @param x coordinates where D is computed
/// @param D constitutive matrix
void TPZElasticity2D::ComputeD(const TPZVec<REAL> &x, TPZFMatrix<STATE> &D) {
    if (fElasticity) {
        TPZManVector<STATE,2> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(x, result, Dres);
        REAL E = result[0];
        REAL nu = result[1];
        TPZLinearElasticityConstitutive loc(E, nu);
        loc.SetPlaneStress(fPlaneStress);
        loc.SetDimension(2);
        loc.ComputeStiffnessMatrix(D);
    } else {
        //D matrix already computed with fE_def and fnu_def
        fConstitutiveLaw.ComputeStiffnessMatrix(D);
    }
}

void TPZElasticity2D::Contribute(const TPZMaterialDataT<STATE> &data,
                                 REAL weight,
                                 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShape(data,weight,ek, ef);
        return;
    }
    if(fUseVecShape == true){
        TPZMaterialDataT<STATE> modifieddata = data;
        // Convert scalar shape functions into vectorial shape functions
        auto &phi = modifieddata.fH1.fPhi;
        auto &dphi = modifieddata.dphix;
        const int nshape = phi.Rows();
        phi.Redim(2,2*nshape);
        dphi.Redim(4,2*nshape);
        for(int ishape = 0; ishape < nshape; ishape++) {
            // X direction
            phi(0,2*ishape) = data.fH1.fPhi(ishape,0);
            dphi(0,2*ishape) = data.dphix(0,ishape);
            dphi(1,2*ishape) = data.dphix(1,ishape);
            // Y direction
            phi(1,2*ishape+1) = data.fH1.fPhi(ishape,0);
            dphi(2,2*ishape+1) = data.dphix(0,ishape);
            dphi(3,2*ishape+1) = data.dphix(1,ishape);
        }
        modifieddata.fShapeType = data.EVecShape;
        // Compute contribution with vectorial shape functions
        ContributeVecShape(modifieddata,weight,ek, ef);
        return;
    }
    
	const TPZFMatrix<REAL> &dphi = data.dphix;
	const TPZFMatrix<REAL> &phi = data.fH1.fPhi;
	const TPZFMatrix<REAL> &axes=data.axes;
	
	int64_t phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if(phc != 1 || dphr != 2 || phr != dphc){
		PZError << "\nTPZElasticity2D.contr, inconsistent input data : \n" <<
		"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
		"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
		return;
		//		PZError.show();
	}


#ifdef PZDEBUG
    if(fForcingFunction && fAnisotropicForcingFunction) {
        PZError << "TPZElasticity2D::Contribute error both forcing functions are set\n";
        DebugStop();
    }
#endif
    TPZFNMatrix<9,STATE> D(3,3);
    ComputeD(data.x, D);
    TPZManVector<STATE,3> floc(ff);
    if(fAnisotropicForcingFunction) {
        TPZManVector<STATE,3> res(3,0.);
        fAnisotropicForcingFunction(data.x, fConstitutiveLaw, res);
        floc[0] = res[0];
        floc[1] = res[1];
    }
	else if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE,3> res(3,0.);
		fForcingFunction(data.x,res);
		floc[0] = res[0];
		floc[1] = res[1];
		floc[2] = res[2];
	}
	
    

    
    TPZFNMatrix<60,STATE> dphidx(3,phr),BMat(3,2*phr,0.),DBMat(3,2*phr);
    TPZAxesTools<STATE>::Axes2XYZ(dphi,dphidx,axes);

    for(int i=0; i<phr; i++) {
        BMat(0,2*i) = dphidx(0,i);
        BMat(2,2*i) = dphidx(1,i);
        BMat(1,2*i+1) = dphidx(1,i);
        BMat(2,2*i+1) = dphidx(0,i);
    }
    D.Multiply(BMat,DBMat);
    ek.AddContribution(0,0,BMat,1,DBMat,0,weight);

	
    for( int in = 0; in < phr; in++ ) {
        for (int col = 0; col < efc; col++)
        {
            ef(2*in, col) += weight * (floc[0]*phi(in,0) - dphidx(0,in)*fPreStressXX - dphidx(1,in)*fPreStressXY);  // direcao x
            ef(2*in+1, col) += weight * (floc[1]*phi(in,0) - dphidx(0,in)*fPreStressXY - dphidx(1,in)*fPreStressYY);// direcao y <<<----
        }
    }
}

void TPZElasticity2D::ContributeBC(const TPZMaterialDataT<STATE> &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,
                                   TPZBndCondT<STATE> &bc) {
    
    
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShapeBC(data,weight,ek, ef,bc);
        return;
    }
    if(fUseVecShape == true){
        TPZMaterialDataT<STATE> modifieddata = data;
        // Convert scalar shape functions into vectorial shape functions
        auto &phi = modifieddata.phi;
        auto &dphi = modifieddata.dphix;
        const int nshape = phi.Rows();
        phi.Redim(2,2*nshape);
        dphi.Redim(4,2*nshape);
        for(int ishape = 0; ishape < nshape; ishape++) {
            // X direction
            phi(0,2*ishape) = data.phi(ishape,0);
            // dphi(0,2*ishape) = data.dphix(0,ishape);
            // dphi(1,2*ishape) = data.dphix(1,ishape);
            // Y direction
            phi(1,2*ishape+1) = data.fH1.fPhi(ishape,0);
            // dphi(2,2*ishape+1) = data.dphix(0,ishape);
            // dphi(3,2*ishape+1) = data.dphix(1,ishape);
        }
        modifieddata.fShapeType = data.EVecShape;
        // Compute contribution with vectorial shape functions
        ContributeVecShapeBC(modifieddata,weight,ek, ef, bc);
        return;
    }

	const TPZFMatrix<REAL> &phi = data.phi;
     int dim = Dimension();

	const auto &BIGNUMBER  = TPZMaterial::fBigNumber;
    
	int phr = phi.Rows();
	short in,jn;
    auto  bcLoadCases = dynamic_cast<TPZMatLoadCasesBC<STATE>&>(bc);
    if (ef.Cols() != bcLoadCases.NumLoadCases()) {
        DebugStop();
    }
	
//		In general when the problem is  needed to stablish any convention for ContributeBC implementations

    int nstate = NStateVariables();

    const auto nloads = this->fNumLoadCases;
    constexpr int nvars = 2;
    const auto &bcNumLoads =
        dynamic_cast<TPZMatLoadCasesBC<STATE>&>(bc);

    TPZFNMatrix<30,STATE> v1(nvars,nvars);
    v1 = bc.Val1();
    TPZManVector<STATE,10> v2(nvars*nloads), dirdisp(2,0.);
    if(bc.HasForcingFunctionBC() && nloads != 1) DebugStop();
    if(bc.HasForcingFunctionBC() && (bc.Type() == 1 || bc.Type() == 2 || bc.Type() == 3)) {
        TPZManVector<STATE> u_ex(2,0.);
        TPZFNMatrix<4,STATE> dsol(2,2,0.);
        bc.ForcingFunctionBC()(data.x,u_ex,dsol);
        TPZFNMatrix<9,STATE> D(3,3);
        ComputeD(data.x, D);
        TPZFNMatrix<3,STATE> strain(3,1,0.);
        fConstitutiveLaw.TensorToVoigt(dsol,strain,true);
        TPZFNMatrix<3,STATE> stress(3,1,0.);
        D.Multiply(strain,stress);
        // Now set the prestress values
        stress(0,0) += fPreStressXX;
        stress(1,0) += fPreStressYY;
        stress(2,0) += fPreStressXY;
        TPZManVector<REAL,3> normal(3,0.);
        normal[0] = data.axes(0,1);
        normal[1] = -data.axes(0,0);
        normal[2] = 0.;
        STATE traction_x = stress(0,0)*normal[0] + stress(2,0)*normal[1];
        STATE traction_y = stress(2,0)*normal[0] + stress(1,0)*normal[1];
        if(bc.Type() == 1) {
            // Neumann or Mixed BC
            v2[0] = traction_x;
            v2[1] = traction_y;      
        } else if(bc.Type() == 2) {
            // Mixed BC
            v2[0] = traction_x + v1(0,0)*u_ex[0] + v1(0,1)*u_ex[1];
            v2[1] = traction_y + v1(1,0)*u_ex[0] + v1(1,1)*u_ex[1];

        } else if(bc.Type() == 3) {
            // Directional Dirichlet BC
            STATE u_dir = bc.Val2()[0]*u_ex[0] + bc.Val2()[1]*u_ex[1];
            v2[0] = BIGNUMBER * u_dir*bc.Val2()[0]+ traction_x;
            v2[1] = BIGNUMBER * u_dir*bc.Val2()[1] + traction_y;
        }
        //std::cout << "traction_x = " << traction_x << " traction_y = " << traction_y << std::endl;
    } else if(bc.HasForcingFunctionBC() && bc.Type() == 0) {
        TPZManVector<STATE> u_ex(2,0.);
        TPZFNMatrix<4,STATE> dsol(2,2,0.);
        bc.ForcingFunctionBC()(data.x,u_ex,dsol);
        v2 = u_ex;
    }
	[&bc = std::as_const(bc),
     &bcNumLoads = std::as_const(bcNumLoads),
     &data = std::as_const(data),
     nvars,nloads]( TPZFMatrix<STATE> &v1, TPZVec<STATE> &v2) {
        if(bc.HasForcingFunctionBC()){
            // v2 is already set
            // v1 = bc.Val1();
        }else {
            for(auto l = 0; l < nloads; l++){
                const auto &val2 = bcNumLoads.GetBCRhsVal(l);
                for(auto i = 0; i < nvars; i++)
                    v2[nvars*l+i] = val2[i];
            }
            v1 = bc.Val1();
        }
    }(v1,v2);

    switch (bc.Type()) {
        case 0 :			// Dirichlet condition
        {
            for(in = 0 ; in < phr; in++) {
                for (int il = 0; il<NumLoadCases(); il++)
                {
                    ef(2*in,il)   += BIGNUMBER * v2[2*il+0] * phi(in,0) * weight;        // forced v2 displacement
                    ef(2*in+1,il) += BIGNUMBER * v2[2*il+1] * phi(in,0) * weight;        // forced v2 displacement
                }
                for (jn = 0 ; jn < phi.Rows(); jn++)
                {
                    ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) *phi(jn,0) * weight;
                }
            }
        }
            break;
            
        case 1 :		// Neumann condition
        {
            for (in = 0; in < phr; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    ef(2*in,il) += v2[2*il+0] * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) +=  v2[2*il+1] * phi(in,0) * weight;      // force in y direction
                }
            }
        }
            break;
            
        case 2 :		// Mixed Condition
        {

            for(in = 0 ; in < phi.Rows(); in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
//                    const auto &v2 = bcLoadCases.GetBCRhsVal(il);
                    ef(2*in,il) += v2[2*il+0] * phi(in,0) * weight;        // force in x direction
                    ef(2*in+1,il) += v2[2*il+1] * phi(in,0) * weight;      // forced in y direction
                }
                
                for (jn = 0 ; jn < phi.Rows(); jn++) {
                    ek(2*in,2*jn) += v1(0,0) * phi(in,0) * phi(jn,0) * weight;         // peso de contorno => integral de contorno
                    ek(2*in+1,2*jn) += v1(1,0) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in+1,2*jn+1) += v1(1,1) * phi(in,0) * phi(jn,0) * weight;
                    ek(2*in,2*jn+1) += v1(0,1) * phi(in,0) * phi(jn,0) * weight;
                }
            }   // este caso pode reproduzir o caso 0 quando o deslocamento
            
            break;
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for(in = 0 ; in < phr; in++) {
               ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
               ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                for (jn = 0 ; jn < phr; jn++) {
                    ek(nstate*in+0,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()[0]*bc.Val2()[0];
                    ek(nstate*in+0,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()[0]*bc.Val2()[1];
                    ek(nstate*in+1,nstate*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()[1]*bc.Val2()[1];
                    ek(nstate*in+1,nstate*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * bc.Val2()[1]*bc.Val2()[0];
                }//jn
            }//in
            break;
            
            
        case 4: // stressField Neumann condition
            {
                REAL v2[2];
                for(in = 0; in < dim; in ++)
                {
                    v2[in] =  ( v1(in,0) * data.normal[0] +
                                v1(in,1) * data.normal[1]);
                }
                // The normal vector points towards the neighbour. The negative sign is there to
                // reflect the outward normal vector.
                for(in = 0 ; in < phi.Rows(); in++) {
                    ef(nstate*in+0,0) += v2[0] * phi(in,0) * weight;
                    ef(nstate*in+1,0) += v2[1] * phi(in,0) * weight;
                    //	cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
                    //	cout << "val2:  " << v2[0]  << endl;
                }
            }
            break;
            
        case 5://PRESSAO DEVE SER POSTA NA POSICAO 0 DO VETOR v2
            {
                TPZFNMatrix<2,STATE> res(2,1,0.);
                for(in = 0 ; in < phi.Rows(); in++)
                {
                    for (int il=0; il<NumLoadCases(); il++)
                    {
                        const auto &v2 = bcLoadCases.GetBCRhsVal(il);
                        ef(nstate*in+0,0) += (v2[0]*data.normal[0]) * phi(in,0) * weight ;
                        ef(nstate*in+1,0) += (v2[0]*data.normal[1]) * phi(in,0) * weight ;
                    }
                    for(jn=0; jn<phi.Rows(); jn++)
                    {
                        for(int idf=0; idf<2; idf++) for(int jdf=0; jdf<2; jdf++)
                        {
                            ek(nstate*in+idf,nstate*jn+jdf) += bc.Val1()(idf,jdf)*data.normal[idf]*data.normal[jdf]*phi(in,0)*phi(jn,0)*weight;
                            //BUG FALTA COLOCAR VAL2
                            //                        DebugStop();
                        }
                    }
                    
                }
            }
            break;
            
            
        }
    } 
}

void TPZElasticity2D::FillDataRequirements(TPZMaterialData &data) const 
{
    data.fNeedsSol = false;
    data.fNeedsNormal = false;

}

void TPZElasticity2D::FillBoundaryConditionDataRequirements(int type,
                                                           TPZMaterialData &data) const
{
    data.fNeedsSol = false;
    data.fNeedsNormal = false;
    if (type == 4 || type == 5 || type == 6) {
        data.fNeedsNormal = true;
    }
}

void TPZElasticity2D::VariableNames(TPZVec<std::string> &names) const {
    names.Resize(31);
    names[0] = "displacement";

	names[1] = "Displacement";
	names[2] = "DisplacementMem";
	names[3] = "Pressure";
	names[10] = "MaxStress";
	names[16] = "PrincipalStress1";
	names[6] = "PrincipalStress2";
	names[7] = "SigmaX";
	names[8] = "SigmaY";
	names[9] = "TauXY";
	names[4] = "Strain";
	names[11] = "SigmaZ";
    names[12] = "sig_x";
    names[13] = "sig_y";
    names[14] = "tau_xy";
    names[15] = "Displacement6";
    names[5] = "Stress";
    names[18] = "J2";
    names[19] = "I1";
    names[20] = "J2Stress";
    names[21] = "I1Stress";
    names[22] = "Alpha";
    names[23] = "PlasticSqJ2";
    names[24] = "PlasticSqJ2El";
    names[25] = "YieldSurface";
    names[26] = "NormalStress";
    names[27] = "ShearStress";
    names[28] = "NormalStrain";
    names[29] = "ShearStrain";
    names[30] = "Young_Modulus";
    names[17] = "Poisson";
}

int TPZElasticity2D::VariableIndex(const std::string &name) const
{
        
	if(!strcmp("displacement",name.c_str()))     return 9;
	if(!strcmp("Displacement",name.c_str()))     return 9;
	if(!strcmp("DisplacementMem",name.c_str()))     return 9;
	if(!strcmp("Pressure",name.c_str()))         return 1;
	if(!strcmp("MaxStress",name.c_str()))        return 2;
	if(!strcmp("PrincipalStress1",name.c_str())) return 3;
	if(!strcmp("PrincipalStress2",name.c_str())) return 4;
	if(!strcmp("SigmaX",name.c_str()))           return 5;
	if(!strcmp("SigmaY",name.c_str()))           return 6;
	if(!strcmp("TauXY",name.c_str()))            return 8;//Cedric
	if(!strcmp("Strain",name.c_str()))           return 11;//Philippe
	if(!strcmp("SigmaZ",name.c_str()))           return 12;//Philippe
    
	if(!strcmp("sig_x",name.c_str()))            return 5;
	if(!strcmp("sig_y",name.c_str()))            return 6;
	if(!strcmp("tau_xy",name.c_str()))           return 8;//Cedric
	if(!strcmp("Displacement6",name.c_str()))    return 7;
	if(!strcmp("Stress",name.c_str()))           return 10;
    if(!strcmp("J2",name.c_str()))           return 20;
    if(!strcmp("I1",name.c_str()))           return 21;
    if(!strcmp("J2Stress",name.c_str()))           return 20;
    if(!strcmp("I1Stress",name.c_str()))           return 21;
    if(!strcmp("Alpha",name.c_str()))        return 22;
    if(!strcmp("PlasticSqJ2",name.c_str()))        return 22;
    if(!strcmp("PlasticSqJ2El",name.c_str()))        return 22;
    if(!strcmp("YieldSurface",name.c_str()))        return 27;
    if(!strcmp("NormalStress",name.c_str()))        return 23;
    if(!strcmp("ShearStress",name.c_str()))        return 24;
    if(!strcmp("NormalStrain",name.c_str()))        return 25;
    if(!strcmp("ShearStrain",name.c_str()))        return 26;
    if(!strcmp("Young_Modulus",name.c_str()))        return 28;
    if(!strcmp("Poisson",name.c_str()))        return 29;
    
    
    
	//   cout << "TPZElasticityMaterial::VariableIndex Error\n";
	return TBase::VariableIndex(name);
}

int TPZElasticity2D::NSolutionVariables(int var) const
{

	switch(var) {
		case 0:
			return 2;
		case 1:
		case 2:
			return 1;
		case 3:
		case 4:
			return 2;
		case 5:
		case 6:
		case 8:
			return 1;
		case 7:
			return 6;
		case 9:
			return 3;
		case 10 : //Stress Tensor
			return 3;
        case 11 : //Strain Tensor
            return 3;
            // SigZ
        case 12:
            return 1;
        case 20:
            return 1;
        case 21:
            return 1;
        case 22:
            return 1;
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
            return 3;
        case 28:
        case 29:
            return 1;
		default:
			return TPZMaterial::NSolutionVariables(var);
	}  
}


void TPZElasticity2D::Solution(const TPZMaterialDataT<STATE> &data,
                               int var, TPZVec<STATE> &Solout)
{
    int numbersol = data.dsol.size();
    int ipos = 0;
    if (fPostProcIndex < numbersol) {
        ipos = fPostProcIndex;
    }
    
    REAL E(fE_def), nu(fnu_def);

    TPZFNMatrix<9,STATE> D(3,3);
    
    if (fElasticity) {
        TPZManVector<STATE, 2> result(2);
        TPZFNMatrix<4, STATE> Dres(0, 0);
        fElasticity(data.x, result, Dres);
        E = result[0];
        nu = result[1];
#ifdef PZDEBUG
        if(fConstitutiveLaw.IsIsotropic() == false)
        {
            DebugStop();
        }
#endif
        TPZLinearElasticityConstitutive loc(E, nu);
        loc.SetPlaneStress(fPlaneStress);
        loc.SetDimension(2);
        loc.ComputeStiffnessMatrix(D);
    } else {
        //D matrix already computed with fE_def and fnu_def
        fConstitutiveLaw.ComputeStiffnessMatrix(D);
    }

    if(var == 28)
    {
        Solout[0] = E;
        return ;
    }
    if(var == 29)
    {
        Solout[0] = nu;
        return;
    }

    
    REAL Eover1MinNu2 = E/(1-nu*nu);
    REAL Eover21PlusNu = E/(2.*(1+nu));
    

    const TPZVec<STATE> &Sol = data.sol[ipos];
    const TPZFMatrix<STATE> &DSol = data.dsol[ipos];
    const TPZFMatrix<REAL> &axes = data.axes;
    TPZFNMatrix<4,STATE> DSolxy(3,2);
    TPZAxesTools<STATE>::Axes2XYZ(DSol,DSolxy,axes);
	TPZFNMatrix<3,STATE> voigtstrain(3,1), voigtstress(3,1);
    voigtstrain(0,0) = DSolxy(0,0); //epsx
    voigtstrain(1,0) = DSolxy(1,1); //epsy
    voigtstrain(2,0) = (DSolxy(1,0)+DSolxy(0,1)); //gamma<sub>xy</sub>
    D.Multiply(voigtstrain,voigtstress);
	REAL epsx = DSolxy(0,0);// du/dx
	REAL epsy = DSolxy(1,1);// dv/dy
	REAL epsxy = 0.5*(DSolxy(1,0)+DSolxy(0,1));
    REAL epsz = fConstitutiveLaw.ComputeStrainZ(voigtstrain);
	REAL SigX = voigtstress(0,0) + fPreStressXX;
	REAL SigY = voigtstress(1,0) + fPreStressYY;
    REAL SigZ = fConstitutiveLaw.ComputeSigmaZ(voigtstrain) + fPreStressZZ;
	REAL TauXY = voigtstress(2,0) + fPreStressXY;
    REAL Pressure = (SigX + SigY + SigZ)/3.;
    REAL aux,Sig1,Sig2,angle;
    
    
    epsx = DSolxy(0,0);// du/dx
    epsy = DSolxy(1,1);// dv/dy
    epsxy = 0.5*(DSolxy(1,0)+DSolxy(0,1));
    
    // REAL lambda = GetLambda(E,nu);
    // REAL mu = GetMU(E,nu);


    
	switch(var) {
		case 0:
			//numvar = 2;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			break;
		case 7:
			//numvar = 6;
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			Solout[3] = 0.;
			Solout[4] = 0.;
			Solout[5] = 0.;
			break;
		case 1:
		case 2:
		case 3:
		case 4:
		case 5:
		case 6:
		case 8:
		case 10:
			
			//numvar = 1;
			Solout[0] = SigX+SigY+SigZ;
            // Pressure variable
			if(var == 1) {
				Solout[0] = (SigX+SigY+SigZ)/3.;
				return;
			}
            // TauXY variable
			if(var == 8) {
				Solout[0] = TauXY;
				return;
			}
			if(var ==5) {
				Solout[0] = SigX;
				return;
			}
			if(var == 6) {
				Solout[0] = SigY;
				return;
			}
			aux = sqrt(0.25*(SigX-SigY)*(SigX-SigY)
					   +(TauXY)*(TauXY));
			// Philippe 13/5/99
			//         if(abs(Tau) < 1.e-10 && abs(SigY-SigX) < 1.e-10) angle = 0.;
			if(fabs(TauXY) < 1.e-10 && fabs(SigX-SigY) < 1.e-10) angle = 0.;
			else angle = atan2(2*TauXY,SigX-SigY)/2.;
            // std::cout << "angle (rad) = " << angle << std::endl;
            // std::cout << "2*TauXY = " << 2*TauXY << " SigX - SigY = " << SigX - SigY << std::endl;
			Sig1 = 0.5*(SigX+SigY)+aux;
			Sig2 = 0.5*(SigX+SigY)-aux;
			if(var == 3 ){
				//numvar = 2;
				Solout[0] = Sig1*cos(angle);
				Solout[1] = Sig1*sin(angle);
				return;
			}
			if(var == 4 ) {
				//numvar = 2;
				Solout[0] = -Sig2*sin(angle);
				Solout[1] = Sig2*cos(angle);
				return;
			}
			if(var == 2) {
				REAL sigmax;
				sigmax = (fabs(Sig1) < fabs(Sig2))? fabs(Sig2) : fabs(Sig1);
				Solout[0] = sigmax;
				return;
			}
			if (var ==10)
			{
				Solout[0] = SigX;
				Solout[1] = SigY;
				Solout[2] = TauXY;
				return;
			}
            std::cout << "Very critical error TPZElasticity2D::Solution\n";
			DebugStop();
			//         Solout[0] /= 0.;
			break;
		case 9:
			Solout[0] = Sol[0];
			Solout[1] = Sol[1];
			Solout[2] = 0.;
			break;
        case 11:
            Solout[0] = epsx;
			Solout[1] = epsy;
			Solout[2] = 2.*epsxy;
            break;
        case 12:
            Solout[0] = SigZ;
            break;
            
        case 20:
        {
            REAL J2 = (SigX-Pressure)*(SigX-Pressure)+(SigY-Pressure)*(SigY-Pressure)+(SigZ-Pressure)*(SigZ-Pressure)+2*TauXY*TauXY;
            J2 /= 2.;
        //    REAL J2 = (pow(SigX + SigY,2) - (3*(-pow(SigX,2) - pow(SigY,2) + pow(SigX + SigY,2) - 2*pow(TauXY,2)))/2.)/2.;
            
            Solout[0]=J2;
            break;
        }
        case 21:
        {
            REAL I1 = SigX+SigY;
            Solout[0]=I1;
            break;
        }
        case 22:
            Solout[0] = 0.;
            break;
        case 23:
            // normal stress
            Solout[0] = SigX;
            Solout[1] = SigY;
            Solout[2] = SigZ;
            break;
        case 24:
            // shear stress
            Solout[0] = TauXY;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 25:
            Solout[0] = epsx;
            Solout[1] = epsy;
            Solout[2] = epsz;
            break;
        case 26:
            Solout[0] = epsxy;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
        case 27:
            Solout[0] = 0.;
            Solout[1] = 0.;
            Solout[2] = 0.;
            break;
		default:
			TBase::Solution(data,var,Solout);
			break;
	}
}


TPZMaterial* TPZElasticity2D::NewMaterial() const
{
    return new TPZElasticity2D();
}


int TPZElasticity2D::ClassId() const
{
    return Hash("TPZElasticity2D") ^
        TBase::ClassId() << 1;
}

void TPZElasticity2D::Read(TPZStream &buf, void *context)
{
    TBase::Read(buf,context);
	buf.Read(&fE_def,1);
	buf.Read(&fnu_def,1);
	buf.Read(&fPreStressXX,1);
	buf.Read(&fPreStressYY,1);
	buf.Read(&fPreStressXY,1);
	buf.Read(&fPreStressZZ,1);
	
	buf.Read(&ff[0],3);
	buf.Read(&fPlaneStress,1);
    buf.Read(&fPostProcIndex);
	
}
	
void TPZElasticity2D::Write(TPZStream &buf, int withclassid) const
{
    TBase::Write(buf,withclassid);
	buf.Write(&fE_def,1);
	buf.Write(&fnu_def,1);
	buf.Write(&fPreStressXX,1);
	buf.Write(&fPreStressYY,1);
	buf.Write(&fPreStressXY,1);
	buf.Write(&fPreStressZZ,1);
	
	buf.Write(&ff[0],3);
	buf.Write(&fPlaneStress,1);
    buf.Write(&fPostProcIndex);

}


void TPZElasticity2D::ContributeVecShape(const TPZMaterialDataT<STATE> &data,
                                         REAL weight,
                                         TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    if(data.fShapeType!=data.EVecShape){
        DebugStop();
    }
    const TPZFMatrix<REAL> &dphi = data.dphix;
	const TPZFMatrix<REAL> &phi = data.fH1.fPhi;
	const TPZFMatrix<REAL> &axes=data.axes;
	
	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
    TPZManVector<STATE,3> floc(ff);
	
    if(fAnisotropicForcingFunction) {
        TPZManVector<STATE,3> res(3,0.);
        fAnisotropicForcingFunction(data.x, fConstitutiveLaw, res);
        floc[0] = res[0];
        floc[1] = res[1];
    }
	else if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE,3> res(3,0.);
		fForcingFunction(data.x,res);
		floc[0] = res[0];
		floc[1] = res[1];
		floc[2] = res[2];
	}

	
    REAL E(fE_def), nu(fnu_def);
    TPZFNMatrix<9,STATE> D(3,3);
    if (fElasticity) {
        TPZManVector<STATE,2> result(2);
        TPZFNMatrix<4,STATE> Dres(0,0);
        fElasticity(data.x, result, Dres);
        E = result[0];
        nu = result[1];
        TPZLinearElasticityConstitutive loc(E, nu);
        loc.SetPlaneStress(fPlaneStress);
        loc.SetDimension(2);
        loc.ComputeStiffnessMatrix(D);
    } else {
        //D matrix already computed with fE_def and fnu_def
        fConstitutiveLaw.ComputeStiffnessMatrix(D);
    }
    

    

	TPZFNMatrix<4,STATE> dphix_i(2,1),dphiy_i(2,1);
	/*
	 * Plain strain materials values
	 */


    TPZFMatrix<REAL> BMat(3,phc,0.),DBMat(3,phc);

	for( int in = 0; in < phc; in++ )
    {
		dphix_i(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);
		dphix_i(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);
		dphiy_i(0,0) = dphi(2,in)*axes(0,0)+dphi(3,in)*axes(1,0);
		dphiy_i(1,0) = dphi(2,in)*axes(0,1)+dphi(3,in)*axes(1,1);

        BMat(0,in) = dphix_i(0,0);
        BMat(1,in) = dphiy_i(1,0);
        BMat(2,in) = dphix_i(1,0) + dphiy_i(0,0);
		
        for (int col = 0; col < efc; col++) 
        {
            ef(in,col) += weight*(   floc[0] * phi(0, in)- dphix_i(0,0)*fPreStressXX - dphix_i(1,0)*fPreStressXY
                                   + floc[1] * phi(1, in)- dphiy_i(1,0)*fPreStressYY - dphiy_i(0,0)*fPreStressXY);
        }		
    }
    // Compute D*B
    DBMat.Zero();
    D.Multiply(BMat, DBMat);
    // Now compute B^T * D * B
    ek.AddContribution(0,0,BMat,1,DBMat,0,weight);
}

void TPZElasticity2D::ContributeVecShapeBC(const TPZMaterialDataT<STATE> &data,
                                           REAL weight,
                                           TPZFMatrix<STATE> &ek,
                                           TPZFMatrix<STATE> &ef,
                                           TPZBndCondT<STATE> &bc) {
    
    const TPZFMatrix<REAL> &phi = data.phi;
    
	const REAL &BIGNUMBER  = TPZMaterial::fBigNumber;
    
	int phc = phi.Cols();
	short in,jn;
        auto  bcLoadCases = dynamic_cast<TPZMatLoadCasesBC<STATE>&>(bc);
    if (ef.Cols() != bcLoadCases.NumLoadCases()) {
        DebugStop();
    }
	
//		In general when the problem is  needed to stablish any convention for ContributeBC implementations

    int nstate = NStateVariables();
    constexpr int nvars = 2;

    const auto nloads = this->fNumLoadCases;
    TPZFNMatrix<4,STATE> v1(nvars,nvars);
    v1 = bc.Val1();
    TPZManVector<STATE,10> v2(nvars*nloads,0.), dirdisp(2,0.);

    const auto &bcNumLoads =
        dynamic_cast<TPZMatLoadCasesBC<STATE>&>(bc);
    if(bcNumLoads.NumLoadCases() != nloads || (nloads != 1 && bc.HasForcingFunctionBC())) {
        DebugStop();
    }
    if(bc.HasForcingFunctionBC() && (bc.Type() == 1 || bc.Type() == 2 || bc.Type() == 3)) {
        TPZManVector<STATE> u_ex(2,0.);
        TPZFNMatrix<4,STATE> dsol(2,2,0.);
        bc.ForcingFunctionBC()(data.x,u_ex,dsol);
        TPZFNMatrix<9,STATE> D(3,3);
        ComputeD(data.x, D);
        TPZFNMatrix<3,STATE> strain(3,1,0.);
        fConstitutiveLaw.TensorToVoigt(dsol,strain,true);
        TPZFNMatrix<3,STATE> stress(3,1,0.);
        D.Multiply(strain,stress);
        // Now set the prestress values
        stress(0,0) += fPreStressXX;
        stress(1,0) += fPreStressYY;
        stress(2,0) += fPreStressXY;
        TPZManVector<REAL,3> normal(3,0.);
        normal[0] = data.axes(0,1);
        normal[1] = -data.axes(0,0);
        normal[2] = 0.;
        STATE traction_x = stress(0,0)*normal[0] + stress(2,0)*normal[1];
        STATE traction_y = stress(2,0)*normal[0] + stress(1,0)*normal[1];
        if(bc.Type() == 1) {
            // Neumann or Mixed BC
            v2[0] = traction_x;
            v2[1] = traction_y;      
        } else if(bc.Type() == 2) {
            // Mixed BC
            v2[0] = traction_x + v1(0,0)*u_ex[0] + v1(0,1)*u_ex[1];
            v2[1] = traction_y + v1(1,0)*u_ex[0] + v1(1,1)*u_ex[1];

        } else if(bc.Type() == 3) {
            // Directional Dirichlet BC
            STATE u_dir = bc.Val2()[0]*u_ex[0] + bc.Val2()[1]*u_ex[1];
            v2[0] = BIGNUMBER * u_dir*bc.Val2()[0]+ traction_x;
            v2[1] = BIGNUMBER * u_dir*bc.Val2()[1] + traction_y;
        }
        //std::cout << "traction_x = " << traction_x << " traction_y = " << traction_y << std::endl;
    } else if(bc.HasForcingFunctionBC() && bc.Type() == 0) {
        TPZManVector<STATE> u_ex(2,0.);
        TPZFNMatrix<4,STATE> dsol(2,2,0.);
        bc.ForcingFunctionBC()(data.x,u_ex,dsol);
        v2 = u_ex;
    }

	[&bc = std::as_const(bc),
     &bcNumLoads = std::as_const(bcNumLoads),
     &data = std::as_const(data),
     nvars,nloads]( TPZFMatrix<STATE> &v1, TPZVec<STATE> &v2) {
        if(bc.HasForcingFunctionBC()){
            // v1 has already been set
        }else {
            for(auto l = 0; l < nloads; l++){
                const auto &val2 = bcNumLoads.GetBCRhsVal(l);
                for(auto i = 0; i < nvars; i++)
                    v2[nvars*l+i] = val2[i];
            }
            v1 = bc.Val1();
        }
    }(v1,v2);

	// auto  bcLoadCases = dynamic_cast<TPZMatLoadCasesBC<STATE>&>(bc);
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phc; in++) {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    ef(in,il) += weight*BIGNUMBER*(v2[2*il+0]*phi(0,in) + v2[2*il+1] * phi(1,in));
                }
				for (jn = 0 ; jn < phc; jn++) {
                    
                    ek(in,jn) += weight*BIGNUMBER*(phi(0,in)*phi(0,jn) + phi(1,in)*phi(1,jn));
				}
			}
			break;
			
		case 1 :			// Neumann condition
            for (in = 0; in < phc; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    ef(in,il)+= weight*(v2[2*il+0]*phi(0,in) + v2[2*il+1]*phi(1,in));
                }
            }
			break;
			
		case 2 :		// condicao mista
			for(in = 0 ; in < phc; in++) 
            {
                for (int il = 0; il <fNumLoadCases; il++) 
                {
                    ef(in,il) += weight * (v2[2*il+0]*phi(0,in) + v2[2*il+1]*phi(1,in));
                }
				
				for (jn = 0; jn <phc; jn++) {
                    
                    ek(in,jn) += v1(0,0)*phi(0,in)*phi(0,jn)*weight 
                    
                                + v1(1,0)*phi(1,in)*phi(0,jn)*weight
                    
                                + v1(0,1)*phi(0,in)*phi(1,jn)*weight
                    
                                + v1(1,1)*phi(1,in)*phi(1,jn)*weight;
				}
			}// este caso pode reproduzir o caso 0 quando o deslocamento
            break;
        case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
            for(in = 0 ; in < phc; in++) {
               ef(in,0) += v2[0] * phi(0,in) * weight;
               ef(in,0) += v2[1] * phi(1,in) * weight;
                for (jn = 0 ; jn < phc; jn++) {
                    for(int i = 0; i < 2; i++) {
                        for(int j = 0; j < 2; j++) {
                            ek(in,jn) += BIGNUMBER * phi(i,in) * phi(j,jn) * weight * bc.Val2()[i]*bc.Val2()[j];
                        }
                    }
                }//jn
            }//in
            break;
    }

}

void TPZElasticity2D::Errors(const TPZMaterialDataT<STATE> &data,
                             TPZVec<REAL> &values) {

    const auto &x = data.x;
    const auto &u = data.sol[0];
    const auto &dudaxes = data.dsol[0];
    const auto &axes = data.axes;
#ifdef PZDEBUG
    if(!this->HasExactSol()){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nThe material has no associated exact solution. Aborting...\n";
        DebugStop();
    }
#endif
    TPZManVector<STATE,3> u_exact(3,0.);
    TPZFNMatrix<9,STATE> du_exact(2,2,0.);
    if(this->HasExactSol()) {
        this->ExactSol()(x,u_exact,du_exact);
    } else {
        DebugStop();
    }
	values[0] = 0.;
	TPZFNMatrix<3,REAL> sigma(3,1,0.),sigma_exact(3,1,0.);
	REAL sigx,sigy,sigxy;
    TPZFNMatrix<6,STATE> dudx(3,dudaxes.Cols());
    TPZAxesTools<STATE>::Axes2XYZ(dudaxes,dudx,axes);
//	du(0,0) = dudx(0,0)*axes(0,0)+dudx(1,0)*axes(1,0);
//	du(1,0) = dudx(0,0)*axes(0,1)+dudx(1,0)*axes(1,1);
//	du(0,1) = dudx(0,1)*axes(0,0)+dudx(1,1)*axes(1,0);
//	du(1,1) = dudx(0,1)*axes(0,1)+dudx(1,1)*axes(1,1);
	
    REAL E(fE_def), nu(fnu_def);

    TPZFNMatrix<9,STATE> D(3,3);
    ComputeD(data.x,D);

    TPZFNMatrix<3,STATE> voigtstrain(3,1), voigtstress(3,1);
    voigtstrain(0,0) = dudx(0,0); //epsx
    voigtstrain(1,0) = dudx(1,1); //epsy
    voigtstrain(2,0) = 0.5*(dudx(1,0)+dudx(0,1)); //gamma<sub>xy</sub>
    D.Multiply(voigtstrain,voigtstress);

	REAL SigX = voigtstress(0,0) + fPreStressXX;
	REAL SigY = voigtstress(1,0) + fPreStressYY;
    REAL SigZ = fConstitutiveLaw.ComputeSigmaZ(voigtstrain) + fPreStressZZ;
	REAL TauXY = voigtstress(2,0) + fPreStressXY;


	//tensoes aproximadas : uma forma
	sigma(0) = SigX;
	sigma(1) = SigY;
	sigma(2) = TauXY;
	
    TPZFNMatrix<3,STATE> eps_exact_voigt(3,1);
    eps_exact_voigt(0,0) = du_exact(0,0); //epsx
    eps_exact_voigt(1,0) = du_exact(1,1); //epsy
    eps_exact_voigt(2,0) = 0.5*(du_exact(1,0)+du_exact(0,1)); //gamma<sub>xy</sub>
    D.Multiply(eps_exact_voigt,sigma_exact);
	//exata
	sigma_exact(0) += fPreStressXX;
	sigma_exact(1) += fPreStressYY;
	sigma_exact(2) += fPreStressXY;
    
	//values[0] = calculo do erro estimado em norma Energia
    values[0] = 0.;
    values[3] = 0.;
    for(int i = 0; i < 3; i++) {
        
        values[0] += (sigma_exact(i)-sigma(i))*(eps_exact_voigt(i)-voigtstrain(i));
        values[3] += sigma_exact(i)*eps_exact_voigt(i);
    }
    if(values[0] < 0 && IsZero(values[0]))
    {
        values[0] = 0.;
    }
	
    if(std::isnan(values[0]) || std::isnan(values[3]) || values[0]<0. || values[3]<0.)
    {
        D.Print("DMat =",std::cout,EMathematicaInput);
        std::cout << "Negative error computed or NaN : " << values[0] << " or " << values[3] << std::endl;
        std::cout << "Sigma approx = " << sigma << std::endl;
        std::cout << "Sigma exact  = " << sigma_exact << std::endl;
        std::cout << "Eps approx = " << voigtstrain << std::endl;
        std::cout << "Eps exact  = " << eps_exact_voigt << std::endl;
        DebugStop();
    }
    
	//values[4] : erro em norma L2 em tensoes
    values[4] = sigma(0)*sigma(0)+sigma(1)*sigma(1)+2*sigma(2)*sigma(2);
    
    //values[5] : erro em norma L2 em sig_xx
    values[5] = sigma(0)*sigma(0);
	
	//values[1] : erro em norma L2 em deslocamentos
	values[1] = (u[0] - u_exact[0])*(u[0] - u_exact[0])+(u[1] - u_exact[1])*(u[1] - u_exact[1]);
	
	values[2] = values[1] + values[0];
}

