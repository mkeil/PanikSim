#include "kernels.h"
#include "deviceFunc.cu"
#include "base.c"
#include <stdio.h>

/* wall force */
__global__ void calcWallForces (float *fwallx, float *fwally, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nw) {
	
	int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
		
	int iw;
	float R = para->R; 
	int InjurySwitch = para->InjurySwitch;
	
	int can_see;
	float tmpr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay ;
	
	if (i <= N) {
        for(iw=0; iw<Nw; iw++) {
			
            WallParticleRelation(iw,i,&tmpr,&can_see,Y[i],X[i],WP,para);
			
            if((can_see==1)&&(tmpr<=R)) {

                /* init */
                tmp_fpsx = tmp_fpsy = 0.0;
                tmp_fyox = tmp_fyoy = 0.0;
                tmp_ftax = tmp_ftay = 0.0;

                /* psychological force */
				
                WallPsychForce(iw,i,tmpr,&tmp_fpsx,&tmp_fpsy, D[i], para);
                /* Young and tangential forces */
                if(tmpr<=0.5*D[i]) {
                    
					WallYoungForce(iw,i,tmpr,&tmp_fyox,&tmp_fyoy, D[i], para);
				
                    WallTangForce_FS1(iw,i,tmpr, &tmp_ftax, &tmp_ftay, D[i], VX[i], VY[i], para);

                }
                /* summing wall forces */
                if(Injured[i]==0) {
                    fwallx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                    fwally[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                } else { /* ie. if Injured[i]=1 */
                    fwallx[i] += tmp_fyox + tmp_ftax;
                    fwally[i] += tmp_fyoy + tmp_ftay;
                }

                /* sum of magnitude of touching forces */
                if(InjurySwitch==1) {
                    ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                }

                // measuring x component of touching force exerted on walls left and right from exit 
				/* lasse ich erstmal weg
                if((iw==1)||(iw==7)) {
                    FW_x -= tmp_fyox + tmp_ftax;
                }
				*/
            }
        }
    }
}

__global__ void calcWPointForces (float *fwpointx, float *fwpointy, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nwp) {
	
	int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
	
	
	float R = para->R; 
	int InjurySwitch = para->InjurySwitch;
	int iwp;
	int can_see;
	float tmpr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay;
	
	if (i <=  N) {
        for(iwp=0; iwp<Nwp; iwp++) {
						
            WPointParticleRelation(iwp,i,&tmpr,&can_see, Y[i], X[i], WP);
            if((can_see==1)&&(tmpr<=R)) {

                /* init */
                tmp_fpsx = tmp_fpsy = 0.0;
                tmp_fyox = tmp_fyoy = 0.0;
                tmp_ftax = tmp_ftay = 0.0;

                /* computing forces */
                WPointPsychForce(iwp,i,tmpr,&tmp_fpsx,&tmp_fpsy, X[i], Y[i], D[i], WP[iwp], para);
                if(tmpr<=0.5*D[i]) {
                    
					WPointYoungForce(iwp,i,tmpr,&tmp_fyox,&tmp_fyoy, X[i], Y[i], D[i], WP[iwp], para);

					WPointTangForce_FS1(iwp,i,tmpr,&tmp_ftax,&tmp_ftay, X[i], Y[i], D[i], VX[i], VY[i], WP[iwp], para);
                }

                /* summing forces */
                if(Injured[i]==0) {
                    fwpointx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                    fwpointy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                } else { /* ie. if Injured[i]=1 */
                    fwpointx[i] += tmp_fyox + tmp_ftax;
                    fwpointy[i] += tmp_fyoy + tmp_ftay;
                }

                /* sum of magnitude of touching forces */
                if(InjurySwitch==1) {
                    ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                }

                // measuring x component of touching force exerted on walls left and right from exit
                // erstmal rausgenommen   
				//
                // if((iwp==0)||(iwp==3)) {
                //    FW_x -= tmp_fyox + tmp_ftax;
                // }

            }
        }
    }
}

__global__ void calcColumnForces (float *fcolx, float *fcoly, float *ftmagsum,float *D, int *Injured, float *X, float *Y, float *VX, float *VY, parameter *para, int N) {
    
	int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
	
	int InjurySwitch = para->InjurySwitch;
	int ColumnSwitch = para->ColumnSwitch;
	float ColumnCenterX = para-> ColumnCenterX;
	float ColumnCenterY = para-> ColumnCenterY;
	float ColumnD = para-> ColumnD;
	
	float A = para-> A;
	float B = para-> B;
	float C_Young = para-> C_Young;
	float Kappa = para-> Kappa;
	
	float R = para -> R;
		
	// lokale Variable
	float tmprsqr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay, rx, ry, f_over_r, scal_prod_over_rsqr, tmpr;
	
	/* 1.4
     * column
     */
	if (i <= N) {
	
		switch(ColumnSwitch) {
		default:
		case 0: {
				// for(i=0; i<N; i++) {
					fcolx[i] = fcoly[i] = 0.0;
				// }
				 break;
			}
		case 1: {
				// for(i=0; i<N; i++) {
					tmprsqr = SQR(X[i]-ColumnCenterX)+SQR(Y[i]-ColumnCenterY);
					if(tmprsqr<=SQR(R)) {
						tmpr=sqrt(tmprsqr);

						/* init */
						tmp_fpsx = tmp_fpsy = 0.0;
						tmp_fyox = tmp_fyoy = 0.0;
						tmp_ftax = tmp_ftay = 0.0;

						/* computing forces */
						/* psychological */
						f_over_r = A * exp(-(tmpr-0.5*(D[i]+ColumnD))/B) / tmpr;
						tmp_fpsx = (X[i]-ColumnCenterX) * f_over_r;
						tmp_fpsy = (Y[i]-ColumnCenterY) * f_over_r;
						/* touching */
						if(tmpr<=0.5*(D[i]+ColumnD)) {
							/* Young */
							f_over_r = 2.0*C_Young*(0.5*(D[i]+ColumnD)-tmpr) / tmpr;
							tmp_fyox = (X[i]-ColumnCenterX) * f_over_r;
							tmp_fyoy = (Y[i]-ColumnCenterY) * f_over_r;
							/* friction */
							rx = X[i]-ColumnCenterX;
							ry = Y[i]-ColumnCenterY;
							scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(tmpr);

							tmp_ftax =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr)
										 * (   ry * scal_prod_over_rsqr );
							tmp_ftay =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr)
										 * ( - rx * scal_prod_over_rsqr );


						}


						/* summing forces */
						if(Injured[i]==0) {
							fcolx[i] = tmp_fpsx + tmp_fyox + tmp_ftax;
							fcoly[i] = tmp_fpsy + tmp_fyoy + tmp_ftay;
						} else { /* ie. if Injured[i]==1 */
							fcolx[i] = tmp_fyox + tmp_ftax;
							fcoly[i] = tmp_fyoy + tmp_ftay;
						}


						/* sum of magnitude of touching forces */
						if(InjurySwitch==1) {
							ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
						}
					}
				// }
				break;
			}
		}
	}
}

__global__ void calcInjuryForces (float *fsmokex, float *fsmokey, float *VX, float *VY, float *V0of, int *Injured, float *ftmagsum, int N, int UpdNum, float *SimTime, float *Phi, float *X, float *D, parameter *para){

	// lokale Variablen
	float x_smokefront, tmpf,tmpr;
	
	int InjurySwitch = para -> InjurySwitch; 
	float FCrush_over_1m = para -> FCrush_over_1m;
	float SmokeStartTime = para -> SmokeStartTime;
	float VSmoke = para -> VSmoke;
	float R = para -> R ;
	float A_fire = para -> A_fire;
	float B_fire = para -> B_fire;
	
	int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
	
	
	
	if (i <= N) {
	
		switch(InjurySwitch) {
		case 0: {
				break;
			}
		case 1: {

				/* case: people crushed */
			
				// frisch verletzt
				if((ftmagsum[i]>FCrush_over_1m*PI*D[i])&&(Injured[i]==0)) {
					Injured[i] = 1;
					// NInjured++; wird anschließend neu berechnet
					V0of[i] = 0.0;
				}
				
				break;
			}
		case 2:
		case 3: {

				/* case: smoke front */
				if(SimTime[UpdNum]>=SmokeStartTime) {
					x_smokefront = (SimTime[UpdNum]-SmokeStartTime)*VSmoke;

					
					/* checking position compared to smoke front */
					tmpr = X[i] - x_smokefront;

					/* center of particle behind smoke front: injured */
					
					printf("Index: %d verletzt: %d \n", i, Injured[i]);
					if( tmpr < 0.5*D[i] ) {
						if(Injured[i]==0) {
							Injured[i] = 1;
							
							V0of[i] = 0.0;
							VX[i] = VY[i] = 0.0;
						}
					}
					/* ahead of front but within its interaction range:
					trying to escape */
					if( (tmpr>=0.5*D[i])&&(tmpr<=R) ) {
						tmpf = A_fire*exp(-(tmpr-0.5*D[i])/B_fire);
						fsmokex[i] += cos(Phi[i])*tmpf;
						fsmokey[i] += sin(Phi[i])*tmpf;
					}
					
				}
				break;
			}
		}
	}
}

__global__ void storeOldValues (float* Xprev_d,float* X_d, float* Yprev_d, float*Y_d, int N)
{

    int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
	
	
    if (i <= N) {
        Xprev_d[i] = X_d[i];
        Yprev_d[i] = Y_d[i];
    }
}

__global__ void calcNewValues (float* X_d,float* Y_d,float* VY_d,float* VX_d, float tstep, int N)
{
    int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
	
	if (i <= N) {
        X_d[i] += VX_d[i] * tstep;
        Y_d[i] += VY_d[i] * tstep;
    }
}

__global__ void sumUp (const int *summanden, const int countElements, int* sum) {
	float summe = 0;
	int i;
	for (i = 0; i < countElements; i++) {
		summe = summe + summanden[i];
	}
	*sum = summe;
}

__global__ void storeNewVelocity (float *VX, float *VY, float *V, float *Vdir,  float *Phi, const float *X, const float *Y, const float *D, wall *W,  const float *vxnew, const float *vynew, parameter *para, const int N, float YS) {

	int b_ID = blockIdx.x; 		   
	int i =  b_ID * blockDim.x + threadIdx.x;
    
	if (i <= N) {
        VX[i] = vxnew[i];
        VY[i] = vynew[i];
        V[i] = sqrt(SQR(VX[i])+SQR(VY[i]));
        Vdir[i] = atan2(VY[i],VX[i]);
        Phi[i] = DirectionOfExit(X[i], Y[i], D[i], YS, para, W);
    }
}