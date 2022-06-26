/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *
*             and soils) using Smoothed Particle Hydrodynamics method              *
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Domain.h"
#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

#include <vector>

#define MIN_PS_FOR_NBSEARCH		1.e-6//TODO: MOVE TO CLASS MEMBER

#include <set>

//https://stackoverflow.com/questions/19240540/dynamically-allocating-array-explain/19240932#19240932
template <typename T>
void Initiate (T ***mat,int row){
	//initialize 10 x 20 array:
	*mat = new T*[row];
	for (int i = 0; i < row; i++)
    (*mat)[i] = new T	; //DO NOT FORGET PARENTHESES, PRECEDENCE OF [] is HIGHER THAN DEREF *
	
	// for (int i=0;i<row;i++)
		// cout<< "row size"<<sizeof(*mat[i])/sizeof(float)<<endl;
}

using namespace std;

namespace SPH {
void General(Domain & dom)
{
}


// Constructor
inline Domain::Domain ()
{
    OutputName[0] = "Property1";
    OutputName[1] = "Property2";
    OutputName[2] = "Property3";
    Time    = 0.0;

    Dimension = 2;
    DomSize	= 0.0,0.0,0.0;

    Gravity	= 0.0,0.0,0.0;

    Cellfac = 2.0;

    KernelType	= 0;
    VisEq	= 0;
    Scheme	= 0;
		GradientType = 0;

    XSPH	= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;
    deltatint	= 0.0;
    deltatmin	= 0.0;
    sqrt_h_a = 0.0025;

    TRPR = 0.0;
    BLPF = 0.0;

    GeneralBefore = & General;
    GeneralAfter = & General;


    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
    I = OrthoSys::I;
	
	Vol=0.;
		auto_ts = true;
		
	
	gradKernelCorr = false;


  m_scalar_prop = 0.;
  
  kin_energy_sum = int_energy_sum = 0.;
	
	thermal_solver = false;

  meshcount = 0;
}

inline Domain::~Domain ()
{
	// size_t Max = Particles.Size();
	// for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

inline void Domain::Kernel_Set(Kernels_Type const & KT)
{
	KernelType = KT;
	if (KernelType==2) Cellfac = 3.0; else Cellfac = 2.0;
}

inline void Domain::Viscosity_Eq_Set(Viscosity_Eq_Type const & VQ)
{
	VisEq = VQ;
}

inline void Domain::Gradient_Approach_Set(Gradient_Type const & GT)
{
	GradientType = GT;
}

inline void Domain::AdaptiveTimeStep()
{
	if (deltatint>deltatmin)
	{
		if (deltat<deltatmin)
			deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		else
			deltat		= deltatmin;
	}
	else
	{
		if (deltatint!=deltat)
			deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		else
			deltat		= deltatint;
	}
	
	// if (contact){
		// if (min_force_ts < deltat)
		// //cout << "Step size changed minimum Contact Forcess time: " << 	min_force_ts<<endl;
		// deltat = min_force_ts;
	// }

	// if (deltat<(deltatint/1.0e5))
		// //cout << "WARNING: Too small time step, please choose a smaller time step initially to make the simulation more stable"<<endl;
		// throw new Fatal("Too small time step, please choose a smaller time step initially to make the simulation more stable");
}

inline void Domain::CheckMinTSVel() {
  //Min time step check based on velocity
  // double test	= 0.0;

  // deltatmin	= deltatint;
  // #pragma omp parallel for schedule (static) private(test) num_threads(Nproc)
  // for (int i=0; i<Particles.Size(); i++) {
    // if (Particles[i]->IsFree) {
      // test = 0.4 * Particles[i]->h/(Particles[i]->Cs + norm(Particles[i]->v));
      // if (deltatmin > test ) {
        // omp_set_lock(&dom_lock);
          // deltatmin = test;
        // omp_unset_lock(&dom_lock);
      // }
    // }
  // }
}

int calcHalfPartCount(const double &r, const double &R, const int xinc){
	int ypartcount = -1;
	if ( xinc > 0 ){
		ypartcount = 1;
		double yp = r;
		double xp = r + (double)(xinc - 1 ) *2.*r; 
		double rad = sqrt(yp*yp + xp*xp);
		while( rad <= R -r ){
			yp += 2.*r;
			rad = sqrt(yp*yp + xp*xp);
			ypartcount++;
		}
		ypartcount-=1;
	}
	return ypartcount;
}

inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									double r, double Density, double h, bool Fixed, bool ghost) { //ghost refers to symmetry at bottom z coordinate

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    //size_t PrePS = Particles.Size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy=1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	
	//// GHOST THING

	int ghost_rows = 3; 

	int xy_ghost_part_count[ghost_rows];
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	
	int id_part=0;
  
  std::vector <Vec3_t> x;
  std::vector <double> dens;
	
  if (Dimension==3) {
		int part_per_row = 0;
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2)/*+r*/;
		//Calculate row count for non ghost particles
		while (zp <= (V(2)+Lz -r)){
			k++; 
      zp += 2.*r;      
		}
		cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
    //Allocate
		
		k = 0;zp = V(2)/*+r*/;

		while (zp <= (V(2)+Lz -r)) {
			j = 0;
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else    
					//Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));
          x.push_back(Vec3_t(xp,yp,zp));
          dens.push_back(Density);
					if (zp == V(2))
						part_per_row++;
					id_part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
      zp += 2.*r;
		}
		cout << "Particles per row: "<<part_per_row<<endl;
    
    m_x  = new Vec3_t [x.size()];
	
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /x.size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<x.size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<x.size(); i++)//Like in Domain::Move
		#endif
		{
      m_x[i]    = x[i];
      m_rho[i]  =dens[i]
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;
}

inline void Domain::MoveGhost(){

	for (int gp=0; gp<GhostPairs.Size(); gp++){
		int  i = GhostPairs[gp].first;
		int gi = GhostPairs[gp].second;
		
    //ASSUMING SYMMETRY
		//See normal direction, if it is vertical
    // tg axis is the same speed
		Particles[gi]-> v  = Particles[i]-> v;
    
    int axis = Particles[gi]-> ghost_plane_axis;
    
		Particles[gi]-> v[axis]  = - Particles[i]-> v[axis];


		Particles[gi]-> a = 0.; //TO NOT INFLUENCE TIME STEP
		
	}
}




}; // namespace SPH
