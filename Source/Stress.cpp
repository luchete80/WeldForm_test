
// Not a particular integration scheme, only adding +dt
inline void Domain::CalcStressStrain(double &dt) {
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	StrainRate 	  = FromFlatSym			        (strrate);
	RotationRate 	= FromFlatAntiSymNullDiag	(rotrate);	
  
	// Jaumann rate terms
	Mat3_t RotationRateT, Stress,SRT,RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	double dep = 0.;
  double sig_trial = 0.;

	// if (FirstStep)
		// ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
	// ShearStressb	= ShearStressa;
	// ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;

	ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

	if (Fail == 1) {
		double J2	= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));
		//Scale back, Fraser Eqn 3-53
		sig_trial = sqrt(3.0*J2);
		ShearStress	= std::min((Sigmay/sig_trial),1.0)*ShearStress;
    
    if      (Material_model == HOLLOMON )       Sigmay = mat->CalcYieldStress(pl_strain); 
		else if  (Material_model == JOHNSON_COOK )  Sigmay = mat->CalcYieldStress(pl_strain,eff_strain_rate,T);
      // if (Material_model == BILINEAR ){
        // //Sigmay = Fy0 + pl_strain*Et
      // } else if (Material_model == HOLLOMON ){
        // Sigmay = mat->CalcYieldStress(pl_strain); 
      // }
			
		if ( sig_trial > Sigmay) {
			if (Material_model == HOLLOMON ){
				Et = mat->CalcTangentModulus(pl_strain); //Fraser 3.54
				Et_m = Et;
			}
			else if (Material_model == JOHNSON_COOK ){// //TODO: > BILINEAR				
				Et = mat->CalcTangentModulus(pl_strain, eff_strain_rate, T); //Fraser 3.54
			} else if (Material_model > BILINEAR ) {//Else Ep = 0
        //cout << "Calculating Ep"<<endl;
				Ep = mat->Elastic().E()*Et/(mat->Elastic().E()-Et);
				// if (Ep < 0)
					// cout << "ATTENTION Material Ep <0 "<<Ep<<", Et" << Et <<", platrain"<<pl_strain<<"effstrrate"<<eff_strain_rate<<endl;
			}
			if (Ep<0) Ep = 1.*mat->Elastic().E();
			//Common for both methods
			//if (Ep>0) {
			dep=( sig_trial - Sigmay)/ (3.*G + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			//cout << "dep: "<<dep<<endl;
			pl_strain += dep;
			delta_pl_strain = dep; // For heating work calculation
			//if (Material_model < JOHNSON_COOK ) //In johnson cook there are several fluences per T,eps,strain rate
			if (Material_model == BILINEAR )
			  Sigmay += dep*Ep;
   
			// dep=( sig_trial - Sigmay)/ (3.*G + Ep);	//Fraser, Eq 3-49 TODO: MODIFY FOR TANGENT MODULUS = 0
			// pl_strain += dep;
			// Sigmay += dep*Ep;
    } //plastic
  }//Fail

	//ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);
	Sigma			= -Pressure * OrthoSys::I + ShearStress;	//Fraser, eq 3.32
	
	if ( dep > 0.0 ) {
		double f = dep/Sigmay;
		Strain_pl(0,0) += f*(Sigma(0,0)-0.5*(Sigma(1,1) + Sigma(2,2) ));
		Strain_pl(1,1) += f*(Sigma(1,1)-0.5*(Sigma(0,0) + Sigma(2,2) ));
		Strain_pl(2,2) += f*(Sigma(2,2)-0.5*(Sigma(0,0) + Sigma(1,1) ));
		Strain_pl(0,1) += 1.5*f*(Sigma(0,1));
		Strain_pl(0,2) += 1.5*f*(Sigma(0,2));
		Strain_pl(1,2) += 1.5*f*(Sigma(1,2));
	}	
  
  Strain	= dt*StrainRate + Strain;


	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}