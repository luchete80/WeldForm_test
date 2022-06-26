#include "Domain.h"

///////////////////////////////////////////////////////////////////
/*  NEW FUNCTION TO CALCULATE ACCELERATION, REMOVING STRAIN AND ROTATION RATES
*** AND ALSO DENSITY; IN ORDER TO CALCULATE THEM SEPARATELY *//////
//  NOTE: ONLY FOR FREE PARTICLES
namespace SPH{
inline void Domain::CalcAccel() {
  
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	{
  for (size_t p=0; p<particle_count;p++) {
    
    double h	= ( h[i]+ h[j])/2;
    Vec3_t xij	= P1->x[i] - P2->x[j];

    double rij	= norm(xij);
    
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

    di = rho[i];
    mi = m[i];

    dj = rho[j];
    mj = m[j];

    Vec3_t vij	= v[i] - v[j];
    double GK	= GradKernel(Dimension, KernelType, rij/h, h);
    double K	= Kernel(Dimension, KernelType, rij/h, h);
		
		// double GK	= m_kernel.gradW(rij/h);
		// double K		= m_kernel.W(rij/h);
		
    //m_clock_begin = clock();
		// Artificial Viscosity
    Mat3_t PIij;
    set_to_zero(PIij);
    if (Alpha!=0.0 || Beta!=0.0)
    {
      double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
      double Cij;
      double Ci,Cj;
      Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
      Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
      Cij = 0.5*(Ci+Cj);
      
      if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
    }
    //m_forces_artifvisc_time += (double)(clock() - m_clock_begin) / CLOCKS_PER_SEC;

    Mat3_t Sigmaj,Sigmai;
    set_to_zero(Sigmaj);
    set_to_zero(Sigmai);
    Sigmai = P1->Sigma;
    Sigmaj = P2->Sigma;

    // Tensile Instability
    Mat3_t TIij;
    set_to_zero(TIij);
    if (P1->TI > 0.0 || P2->TI > 0.0) 
      TIij = pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);
      //TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);
		

		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		

		// NEW
		if (GradientType == 0)
			Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
		else
			Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp);
		} else {


				P1->a					+= mj * temp;

				P2->a					-= mi * temp;

  }//MAIN FOR IN PART

}

//Similar but not densities
inline void Domain::CalcRateTensors() {
  Particle *P1, *P2;
  
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
  for (size_t i=0; i<SMPairs[k].Size();i++) {
    P1	= Particles[SMPairs[k][i].first];
    P2	= Particles[SMPairs[k][i].second];	
    
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	//Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

  double clock_begin;

	// if ((rij/h)<=Cellfac)
	// {
  double di=0.0,dj=0.0,mi=0.0,mj=0.0;

    di = P1->Density;
    mi = P1->Mass;

    dj = P2->Density;
    mj = P2->Mass;

		Vec3_t vij	= P1->v - P2->v;
		
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);
	

		Mat3_t Sigmaj,Sigmai;
		set_to_zero(Sigmaj);
		set_to_zero(Sigmai);
		Sigmai = P1->Sigma;
		Sigmaj = P2->Sigma;

		// NoSlip BC velocity correction
		Vec3_t vab = vij;

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);
		
		//NEW
		Mat3_t GKc[2];
		double m, mc[2];
		GKc[0] = GK * P1->gradCorrM;
		GKc[1] = GK * P2->gradCorrM;
		if (gradKernelCorr){
		}

    //m_clock_begin = clock();		
		
    Mat3_t StrainRate_c[2],RotationRate_c[2]; //Corrected gradients

		// // // // Calculation strain rate tensor
		// // // //ORIGINAL FORM			
		//if (!gradKernelCorr){
    StrainRate(0,0) = 2.0*vab(0)*xij(0);
    StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
    StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
    StrainRate(1,0) = StrainRate(0,1);
    StrainRate(1,1) = 2.0*vab(1)*xij(1);
    StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
    StrainRate(2,0) = StrainRate(0,2);
    StrainRate(2,1) = StrainRate(1,2);
    StrainRate(2,2) = 2.0*vab(2)*xij(2);
    StrainRate	= -0.5 * GK * StrainRate;
    
    RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
    RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
    RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
    RotationRate(1,0) = -RotationRate(0,1);
    RotationRate(2,0) = -RotationRate(0,2);
    RotationRate(2,1) = -RotationRate(1,2);
    RotationRate	  = -0.5 * GK * RotationRate;
    
			// if (StrainRate(2,2)<-1.e-3)

			Mat3_t gradv[2],gradvT[2];
			
			//cout<<"gradv"<<gradv[0]<<endl;
			Vec3_t gradK; 
			Mult(GK * P1->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[0]); //outer product. L, velocity gradient tensor
			Mult(GK * P2->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[1]); //outer product. L, velocity gradient tensor
			
			for (int i=0;i<2;i++){
				Trans(gradv[i],gradvT[i]);
				StrainRate_c[i] 	= -0.5*(gradv[i] + gradvT[i]);
				RotationRate_c[i] = -0.5*(gradv[i] - gradvT[i]);
			}
      
		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		
		if (gradKernelCorr){
			for (int i=0;i<2;i++){
				Mult (GKc[i], xij, vc[i]);
			}
		}


		if (Dimension == 2) temp(2) = 0.0;
		
		if (!gradKernelCorr){
			temp1 = dot( vij , GK*xij );
		} else {
			for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
				temp1_c[i] = dot( vij , vc[i] );
			}
		}
    
    clock_begin = clock();
		// Locking the particle 1 for updating the properties
		omp_set_lock(&P1->my_lock);
			// if (!gradKernelCorr){
				// P1->dDensity	+= mj * (di/dj) * temp1;
			// } else{
				// P1->dDensity	+= mj * (di/dj) * temp1_c[0];
			// }		

      float mj_dj= mj/dj;

        P1->StrainRate 		= P1->StrainRate + mj_dj*StrainRate;
        P1->RotationRate 	= P1->RotationRate + mj_dj*RotationRate;

      float mi_di = mi/di;

        P2->StrainRate	 = P2->StrainRate + mi_di*StrainRate;
        P2->RotationRate = P2->RotationRate + mi_di*RotationRate;

    }//FOR PAIRS
  }//FOR NPROC
}

// TODO: USED CALCULATED KERNELKS
inline void Domain::CalcDensInc() {
  Particle *P1, *P2;
  
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
  for (size_t i=0; i<SMPairs[k].Size();i++) {
    P1	= Particles[SMPairs[k][i].first];
    P2	= Particles[SMPairs[k][i].second];	
    
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	//Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

  double clock_begin;

	// if ((rij/h)<=Cellfac)
	// {
  double di=0.0,dj=0.0,mi=0.0,mj=0.0;

    di = P1->Density;
    mi = P1->Mass;

    dj = P2->Density;
    mj = P2->Mass;

		Vec3_t vij	= P1->v - P2->v;
		
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);
	
		//NEW
		Mat3_t GKc[2];
		double m, mc[2];
		GKc[0] = GK * P1->gradCorrM;
		GKc[1] = GK * P2->gradCorrM;
		if (gradKernelCorr){
		}
		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		
		if (gradKernelCorr){
			for (int i=0;i<2;i++){
				Mult (GKc[i], xij, vc[i]);
			}
		}


		if (Dimension == 2) temp(2) = 0.0;
		
		if (!gradKernelCorr){
			temp1 = dot( vij , GK*xij );
		} else {
			for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
				temp1_c[i] = dot( vij , vc[i] );
			}
		}
    
    clock_begin = clock();
		// Locking the particle 1 for updating the properties
		omp_set_lock(&P1->my_lock);
			if (!gradKernelCorr){
				P1->dDensity	+= mj * (di/dj) * temp1;
			} else{
				P1->dDensity	+= mj * (di/dj) * temp1_c[0];
			}		


		omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		omp_set_lock(&P2->my_lock);
			if (!gradKernelCorr){
				P2->dDensity	+= mi * (dj/di) * temp1;							
			}else {
				P2->dDensity	+= mi * (dj/di) * temp1_c[1];
			}
		omp_unset_lock(&P2->my_lock);
    }//FOR PAIRS
  }//FOR NPROC
}


};//SPH