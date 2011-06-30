#pragma once
#include "types.h"

__global__ void calcWallForces (float *fwallx, float *fwally, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nw);
	
__global__ void calcWPointForces (float *fwpointx, float *fwpointy, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nwp); 

__global__ void calcColumnForces (float *fcolx, float *fcoly, float *ftmagsum,float *D, int *Injured, float *X, float *Y, float *VX, float *VY, parameter *para, int N); 

__global__ void calcInjuryForces (float *fsmokex, float *fsmokey, float *VX, float *VY, float *V0of, int *Injured, float *ftmagsum, int N, float SimTime, float *Phi, float *X, float *D, parameter *para);

__global__ void sumForces (float *fsumx,float *fsumy,  float *tStepVector, const float sqrt_fact, const float *VX,const float *VY,const float *V0of,const float *Phi,const float *fpairx,const float *fwallx,const float *fwpointx,const float *fpairy,const float *fwally,const float *fwpointy,const float *fsmokex,const float *fsmokey,const float *fcolx,const float *fcoly, const int N, const parameter *para);

__global__ void NewVelocity (float *vxnew, float *vynew, const float *fsumx, const float *fsumy, const float *VX, const float *VY, const int *Injured, const int N, const float *tStepVector, parameter *para); 

__global__ void storeOldValues (float* Xprev_d,float* X_d, float* Yprev_d, float*Y_d, int N); 

__global__ void calcNewValues (float* X_d,float* Y_d,float* VY_d,float* VX_d, float tstep, int N); 

__global__ void sumUp (const int *summanden, const int countElements, int* sum) ; 

__global__ void getMinTimeStep (const float *tStepVector, const int countElements, float *min);