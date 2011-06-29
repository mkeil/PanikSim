#pragma once

#include <glog/logging.h>
#include <stdio.h>

#include "types.h"
#include "kernels.h"
#include "deviceFunc.h"
#include "update.h"
#include "prepareParameter.h"




void Upd(parameter para_h)
{
	
    int allocN,i,j,k,l,mx,my,m,j_old,j_new;
    float *fwallx,*fwally,*fwpointx,*fwpointy,*fpairx,*fpairy,
          *fspx,*fspy,*fsumx,*fsumy,*vxnew,*vynew,tstep,tmpr,
          tmp_fpsx,tmp_fpsy,tmp_fyox,tmp_fyoy,tmp_ftax,tmp_ftay,
          tmprsqr,sqrt_fact,ksi,eta,vnew,*ftmagsum,*fsmokex,*fsmokey,
          x_smokefront, tmpf, f_over_r, rx, ry, scal_prod_over_rsqr, 
		  *fcolx,*fcoly;



    /* 0 */
    allocN=N;
    fwallx=vector(0,allocN-1);
    fwally=vector(0,allocN-1);
    
	fwpointx=vector(0,allocN-1);
    fwpointy=vector(0,allocN-1);
    
	fpairx=vector(0,allocN-1);
    fpairy=vector(0,allocN-1);
    
	fsmokex=vector(0,allocN-1);
    fsmokey=vector(0,allocN-1);

    fspx=vector(0,allocN-1);
    fspy=vector(0,allocN-1);
    fsumx=vector(0,allocN-1);
    fsumy=vector(0,allocN-1);
    vxnew=vector(0,allocN-1);
    vynew=vector(0,allocN-1);
    ftmagsum=vector(0,allocN-1);
    fcolx=vector(0,allocN-1);
    fcoly=vector(0,allocN-1);

	// printf("alle normalen Pointer wurden alloziert.\n");
	
	
	
	float *Xprev_d, *Yprev_d, *X_d, *Y_d, *VX_d, *VY_d, *fwallx_d, *fwally_d,*fwpointx_d, *fwpointy_d, *ftmagsum_d, *D_d, *fcolx_d, *fcoly_d, *fsmokex_d, *fsmokey_d, *V0of_d, *SimTime_d, *Phi_d;
	int *Injured_d;
	int *NInjured_d; 
	
	wpoint *WP_d;
	
	parameter *para_d;
	 
	
	
	
	
	int sizeFloatVector = N * sizeof(float);

	cudaError_t error; 
    
	error = cudaMalloc (&NInjured_d,sizeof(int)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // Anzahl der Verletzten
	error = cudaMalloc (&Injured_d,N * sizeof(int)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&WP_d,NWP * sizeof(wpoint)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&para_d,sizeof(parameter)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	error = cudaMalloc (&fwpointx_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&fwpointy_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMalloc (&fwallx_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&fwally_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMalloc (&ftmagsum_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&D_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	error = cudaMalloc (&fcolx_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&fcoly_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMalloc (&fsmokex_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMalloc (&fsmokey_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMalloc (&V0of_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&SimTime_d, sizeFloatVector);CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMalloc (&Xprev_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&Yprev_d, sizeFloatVector);CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&X_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&Y_d,sizeFloatVector);CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&VY_d,sizeFloatVector);CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&VX_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

	error = cudaMalloc (&Phi_d,sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	LOG (INFO) << "Device Pointer alloziert.";
	
	
    /* 1 */

    /* 1.0 */
    /* default values */
    tstep = DefaultDeltaT;
    for(i=0; i<N; i++) {
        fwallx[i] = 0.0;
        fwally[i] = 0.0;
        fwpointx[i] = 0.0;
        fwpointy[i] = 0.0;
        fpairx[i] = 0.0;
        fpairy[i] = 0.0;
        fsmokex[i] = 0.0;
        fsmokey[i] = 0.0;

        fspx[i] = 0.0;
        fspy[i] = 0.0;
        fsumx[i] = 0.0;
        fsumy[i] = 0.0;
        ftmagsum[i] = 0.0;
        fcolx[i] = 0.0;
        fcoly[i] = 0.0;
    }
    FW_x=0.0;
	

	// temporäre Kraftvektoren auf dem Device mit 0 initialsieren
	
	
	error = cudaMemset(fwallx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemset(fwally_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMemset(fcolx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemset(fcoly_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

	error = cudaMemset(fsmokex_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemset(fsmokey_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	error = cudaMemset(fwpointx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemset(fwpointy_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	error = cudaMemset(ftmagsum_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	// error = cudaMemcpy(fwallx_d, fwallx, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	// error = cudaMemcpy(fwally_d, fwally, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	// error = cudaMemcpy(fcolx_d, fcolx, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	// error = cudaMemcpy(fcoly_d, fcoly, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	// error = cudaMemcpy(fwpointx_d, fwpointx, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	// error = cudaMemcpy(fwpointy_d, fwpointy, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	// error = cudaMemcpy(ftmagsum_d, ftmagsum, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	LOG (INFO) << "device kraftvektoren initialisiert";
	
	
	// nötige Werte hochkopieren
	
	error = cudaMemcpy(D_d, D, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(Injured_d, Injured, N * sizeof(int), cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

	error = cudaMemcpy(X_d, X, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(Y_d, Y, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
    error = cudaMemcpy(VY_d, VY, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(VX_d, VX, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	error = cudaMemcpy(WP_d, WP, NWP * sizeof(wpoint), cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
		
	error = cudaMemcpy(para_d, &para_h, sizeof(parameter), cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    	
	LOG (INFO) << "Alle Werte wurden hochkopiert";
	
	
	dim3 dimBlock(32); 
	dim3 dimGrid((para_h.N0 + dimBlock.x - 1) / dimBlock.x);
	
	// LOG(WARNING) << "Block dimensions: " << dimBlock.x << " " << dimBlock.y << " " << dimBlock.z;
	// LOG(WARNING) << "Grid dimensions: " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z;
	
	calcWallForces<<<dimGrid, dimBlock>>> (fwallx_d, fwally_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, WP_d, VX_d, VY_d, para_d, N, NW); 
	cudaThreadSynchronize();
	// berechnete Werte calcWallForces zurück
	error = cudaMemcpy(fwallx, fwallx_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(fwally, fwally_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	// ftmagsum wird im nächsten Kernel nochmal verwendet, darum erst danach zurück
	
	calcWPointForces <<<dimGrid, dimBlock>>> (fwpointx_d, fwpointy_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, WP_d, VX_d, VY_d, para_d, N, NWP);
	cudaThreadSynchronize();
	// berechnete Werte calcWPointForces zurück
	error = cudaMemcpy(fwpointx, fwpointx_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(fwpointy, fwpointy_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(ftmagsum, ftmagsum_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


	
    /* 1.1 */
    /* wall force */
    // for(i=0; i<N; i++) {
        // for(iw=0; iw<NW; iw++) {

            // WallParticleRelation(iw,i,&tmpr,&can_see);
            // if((can_see==1)&&(tmpr<=R)) {

                // /* init */
                // tmp_fpsx = tmp_fpsy = 0.0;
                // tmp_fyox = tmp_fyoy = 0.0;
                // tmp_ftax = tmp_ftay = 0.0;

                // /* psychological force */
                // WallPsychForce(iw,i,tmpr,&tmp_fpsx,&tmp_fpsy);
                // /* Young and tangential forces */
                // if(tmpr<=0.5*D[i]) {
                    // WallYoungForce(iw,i,tmpr,&tmp_fyox,&tmp_fyoy);

                    // WallTangForce_FS1(iw,i,tmpr,&tmp_ftax,&tmp_ftay);

                // }
                // /* summing wall forces */
                // if(Injured[i]==0) {
                    // fwallx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                    // fwally[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                // } else { /* ie. if Injured[i]=1 */
                    // fwallx[i] += tmp_fyox + tmp_ftax;
                    // fwally[i] += tmp_fyoy + tmp_ftay;
                // }

                // /* sum of magnitude of touching forces */
                // if(InjurySwitch==1) {
                    // ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                // }

                // /* measuring x component of touching force exerted
                // on walls left and right from exit
                 // -- only in demo mode */

                // if((iw==1)||(iw==7)) {
                    // FW_x -= tmp_fyox + tmp_ftax;
                // }
            // }
        // }
    // }



    /* 1.2 */
    /* wpoint force */
    // for(i=0; i<N; i++) {
        // for(iwp=0; iwp<NWP; iwp++) {

            // WPointParticleRelation(iwp,i,&tmpr,&can_see);
            // if((can_see==1)&&(tmpr<=R)) {

                // /* init */
                // tmp_fpsx = tmp_fpsy = 0.0;
                // tmp_fyox = tmp_fyoy = 0.0;
                // tmp_ftax = tmp_ftay = 0.0;

                // /* computing forces */
                // WPointPsychForce(iwp,i,tmpr,&tmp_fpsx,&tmp_fpsy);
                // if(tmpr<=0.5*D[i]) {
                    // WPointYoungForce(iwp,i,tmpr,&tmp_fyox,&tmp_fyoy);

                    // WPointTangForce_FS1(iwp,i,tmpr,&tmp_ftax,&tmp_ftay);

                // }

                // /* summing forces */
                // if(Injured[i]==0) {
                    // fwpointx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                    // fwpointy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                // } else { /* ie. if Injured[i]=1 */
                    // fwpointx[i] += tmp_fyox + tmp_ftax;
                    // fwpointy[i] += tmp_fyoy + tmp_ftay;
                // }

                // /* sum of magnitude of touching forces */
                // if(InjurySwitch==1) {
                    // ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                // }

                // /* measuring x component of touching force exerted
                   // on walls left and right from exit
                   // -- only in demo mode */

                // if((iwp==0)||(iwp==3)) {
                    // FW_x -= tmp_fyox + tmp_ftax;
                // }

            // }
        // }
    // }



    /* 1.3 */
    /* particle-particle forces */
    for(i=0; i<N; i++) {

        j = (int)floor(X[i]*GX/XS) + G * (int)floor(Y[i]*GY/YS);
        for(k=-1; k<=1; k++) {
            for(l=-1; l<=1; l++) {

                mx = j%G+k;
                my = j/G+l;
                if((mx>=0)&&(mx<GX)&&(my>=0)&&(my<GY)) {

                    m = BIndBd[ (mx+GX)%GX + G * (my%GY) ];
                    /* checking each pair of particles only once */
                    while(m>=i) {
                        m = BInd[m];
                    }
                    if(m!=-1) {
                        do {

                            tmprsqr = SQR(X[i]-X[m]) + SQR(Y[i]-Y[m]);
                            if( tmprsqr <= SQR(R) ) {
                                tmpr = sqrt(tmprsqr);

                                /* init */
                                tmp_fpsx = tmp_fpsy = 0.0;
                                tmp_fyox = tmp_fyoy = 0.0;
                                tmp_ftax = tmp_ftay = 0.0;

                                /* pair forces */
                                /* Force(i,m,...) gives the force exerted by m
                                on i, all forces are symmetric now */
                                PP_PsychForce(i,m,tmpr,&tmp_fpsx,&tmp_fpsy);
                                if(tmpr<=0.5*(D[i]+D[m])) {
                                    PP_YoungForce(i,m,tmpr,&tmp_fyox,&tmp_fyoy);
                                    PP_TangForce_FS1(i,m,tmpr,&tmp_ftax,&tmp_ftay);

                                }

                                /* summing forces */
                                if(Injured[i]==0) {
                                    fpairx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                                    fpairy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                                } else { /* ie. if Injured[i]=1 */
                                    fpairx[i] += tmp_fyox + tmp_ftax;
                                    fpairy[i] += tmp_fyoy + tmp_ftay;
                                }
                                if(Injured[m]==0) {
                                    fpairx[m] -= tmp_fpsx + tmp_fyox + tmp_ftax;
                                    fpairy[m] -= tmp_fpsy + tmp_fyoy + tmp_ftay;
                                } else { /* ie. if Injured[m]=1 */
                                    fpairx[m] -= tmp_fyox + tmp_ftax;
                                    fpairy[m] -= tmp_fyoy + tmp_ftay;
                                }

                                /* sum of magnitude of touching forces */
                                if(InjurySwitch==1) {
                                    ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                                    ftmagsum[m] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                                }
                            }

                            m = BInd[m];
                            while(m>=i) {
                                m = BInd[m];
                            }

                        } while(m!=-1);
                    }
                }
            }
        }
    }
	// benötigte Werte für calcColumnForces hochkopieren
	
	error = cudaMemcpy(ftmagsum_d, ftmagsum, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	calcColumnForces <<<dimGrid, dimBlock>>> (fcolx_d, fcoly_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, VX_d, VY_d, para_d, N); 
	
	CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); error = cudaGetLastError();
	// printf ("Ergebniss von calcColumnForces: %s \n", cudaGetErrorString(error));
	cudaThreadSynchronize();
	
	// berechnete Werte calcColumnForces zurück
	error = cudaMemcpy(fcolx, fcolx_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(fcoly, fcoly_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); 
	
	
	
	
	// benötigte Werte für calcInjuryForces hochkopieren 
	
	error = cudaMemcpy(V0of_d, V0of, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(SimTime_d, SimTime, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(Phi_d, Phi, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
	
	
	calcInjuryForces <<<dimGrid, dimBlock>>> (fsmokex_d, fsmokey_d, VX_d, VY_d, V0of_d, Injured_d,ftmagsum_d, N, UpdNum, SimTime_d, Phi_d, X_d, D_d, para_d);
	cudaThreadSynchronize ();
	
	
	// berechnete Werte calcInjuryForces zurück
	
	error = cudaMemcpy(fsmokex, fsmokex_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(fsmokey, fsmokey_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(VX, VX_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(VY, VY_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(V0of, V0of_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(Injured, Injured_d, N * sizeof(int), cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	error = cudaMemcpy(ftmagsum, ftmagsum_d, sizeFloatVector, cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); 
	
	
	
	// Anzahl der Verletzten neu bestimmen
	// sumUp<<<1,1>>> (Injured_d,N, NInjured_d); cudaThreadSynchronize(); error = cudaGetLastError ();  printf ("Ergebniss sumUp : %s \n",cudaGetErrorString(error));
	
	// error = cudaMemcpy(&NInjured, NInjured_d, sizeof(int), cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
	
    /* 1.4
     * column
     */
    // switch(ColumnSwitch) {
    // default:
    // case 0: {
            // for(i=0; i<N; i++) {
                // fcolx[i] = fcoly[i] = 0.0;
            // }
            // break;
        // }
    // case 1: {
            // for(i=0; i<N; i++) {
                // tmprsqr = SQR(X[i]-ColumnCenterX)+SQR(Y[i]-ColumnCenterY);
                // if(tmprsqr<=SQR(R)) {
                    // tmpr=sqrt(tmprsqr);

                    // /* init */
                    // tmp_fpsx = tmp_fpsy = 0.0;
                    // tmp_fyox = tmp_fyoy = 0.0;
                    // tmp_ftax = tmp_ftay = 0.0;

                    // /* computing forces */
                    // /* psychological */
                    // f_over_r = A * exp(-(tmpr-0.5*(D[i]+ColumnD))/B) / tmpr;
                    // tmp_fpsx = (X[i]-ColumnCenterX) * f_over_r;
                    // tmp_fpsy = (Y[i]-ColumnCenterY) * f_over_r;
                    // /* touching */
                    // if(tmpr<=0.5*(D[i]+ColumnD)) {
                        // /* Young */
                        // f_over_r = 2.0*C_Young*(0.5*(D[i]+ColumnD)-tmpr) / tmpr;
                        // tmp_fyox = (X[i]-ColumnCenterX) * f_over_r;
                        // tmp_fyoy = (Y[i]-ColumnCenterY) * f_over_r;
                        // /* friction */
                        // rx = X[i]-ColumnCenterX;
                        // ry = Y[i]-ColumnCenterY;
                        // scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(tmpr);

                        // tmp_ftax =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr)
                                     // * (   ry * scal_prod_over_rsqr );
                        // tmp_ftay =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr)
                                     // * ( - rx * scal_prod_over_rsqr );


                    // }


                    // /* summing forces */
                    // if(Injured[i]==0) {
                        // fcolx[i] = tmp_fpsx + tmp_fyox + tmp_ftax;
                        // fcoly[i] = tmp_fpsy + tmp_fyoy + tmp_ftay;
                    // } else { /* ie. if Injured[i]==1 */
                        // fcolx[i] = tmp_fyox + tmp_ftax;
                        // fcoly[i] = tmp_fyoy + tmp_ftay;
                    // }


                    // /* sum of magnitude of touching forces */
                    // if(InjurySwitch==1) {
                        // ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                    // }
                // }
            // }
            // break;
        // }
    // }



    /* 1.5 */
    /* injuries */

    // switch(InjurySwitch) {
    // case 0: {
            // break;
        // }
    // case 1: {

            // /* case: people crushed */
            // for(i=0; i<N; i++) {

                // /* newly injured */
                // if((ftmagsum[i]>FCrush_over_1m*PI*D[i])&&(Injured[i]==0)) {
                    // Injured[i] = 1;
                    // NInjured++;
                    // V0of[i] = 0.0;
                // }
            // }
            // break;
        // }
    // case 2:
    // case 3: {

            // /* case: smoke front */
            // if(SimTime[UpdNum]>=SmokeStartTime) {
                // x_smokefront = (SimTime[UpdNum]-SmokeStartTime)*VSmoke;

                // for(i=0; i<N; i++) {
                    // /* checking position compared to smoke front */
                    // tmpr = X[i] - x_smokefront;

                    // /* center of particle behind smoke front: injured */
                    // if( tmpr < 0.5*D[i] ) {
                        // if(Injured[i]==0) {
                            // Injured[i] = 1;
                            // NInjured++;
                            // V0of[i] = 0.0;
                            // VX[i] = VY[i] = 0.0;
                        // }
                    // }
                    // /* ahead of front but within its interaction range:
                    // trying to escape */
                    // if( (tmpr>=0.5*D[i])&&(tmpr<=R) ) {
                        // tmpf = A_fire*exp(-(tmpr-0.5*D[i])/B_fire);
                        // fsmokex[i] += cos(Phi[i])*tmpf;
                        // fsmokey[i] += sin(Phi[i])*tmpf;
                    // }
                // }
            // }
            // break;
        // }
    // }



    /* 2 */

    /* 2.1 preparing update of the eq. of motion */

    sqrt_fact = sqrt(tstep/DefaultDeltaT);
    for(i=0; i<N; i++) {

        /* self-propelling */
        fspx[i] = 1/Tau * (V0of[i]*cos(Phi[i]) - VX[i]);
        fspy[i] = 1/Tau * (V0of[i]*sin(Phi[i]) - VY[i]);

        /* noise */
        if(GaTh!=0.0) {
            ksi = GaussRand(GaMe, GaTh, GaCM);
            eta = 2.0*PI * rand() / (RAND_MAX+1.0);
        } else {
            ksi=0.0;
            eta=0.0;
        }


        /* sum of forces */
        fsumx[i] =   fspx[i] + fpairx[i] + fwallx[i] + fwpointx[i]
                     + sqrt_fact * ksi * cos(eta);
        fsumy[i] =   fspy[i] + fpairy[i] + fwally[i] + fwpointy[i]
                     + sqrt_fact * ksi * sin(eta);


        /* adding smoke force */
        if((InjurySwitch==2)||(InjurySwitch==3)) {
            fsumx[i] += fsmokex[i];
            fsumy[i] += fsmokey[i];
        }
        /* adding force of column */
        switch(ColumnSwitch) {
        default:
        case 0: {
                break;
            }
        case 1: {
                fsumx[i] += fcolx[i];
                fsumy[i] += fcoly[i];
                break;
            }
        }


        /* time step adjustment for velocity change */
        EulTStep( &tstep, sqrt(SQR(fsumx[i])+SQR(fsumy[i])) );


        /* new velocity */
        if(  (Injured[i]==1)
                &&((InjurySwitch==1)||(InjurySwitch==3))
          ) {
            vxnew[i] = 0.0;
            vynew[i] = 0.0;
        } else {
            vxnew[i] = VX[i] + fsumx[i] * tstep;
            vynew[i] = VY[i] + fsumy[i] * tstep;
        }


        /* checking new velocity */
        vnew = sqrt( SQR(vxnew[i]) + SQR(vynew[i]) );
        if(vnew > Vmax) {
            vxnew[i] = vxnew[i]/vnew * Vmax;
            vynew[i] = vynew[i]/vnew * Vmax;
        }
    }




    /* 3 */
    

    // alte Werte hochkopieren
    error = cudaMemcpy(X_d, X, sizeFloatVector, cudaMemcpyHostToDevice); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(Y_d, Y, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(VY_d, VY, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(VX_d, VX, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

	
    storeOldValues <<<dimGrid, dimBlock>>> (Xprev_d,X_d,Yprev_d,Y_d, N);
    cudaThreadSynchronize();

	
    calcNewValues <<<dimGrid, dimBlock>>> (X_d,Y_d,VY_d,VX_d,tstep, N);
    cudaThreadSynchronize();

    // neue Werte zurückkopieren
    error = cudaMemcpy(Xprev, Xprev_d, sizeFloatVector, cudaMemcpyDeviceToHost);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(Yprev, Yprev_d, sizeFloatVector, cudaMemcpyDeviceToHost);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(X, X_d, sizeFloatVector, cudaMemcpyDeviceToHost);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(Y, Y_d, sizeFloatVector, cudaMemcpyDeviceToHost);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    for(i=0; i<N; i++) {

        /* .1 wurde ersetzt durch storeOldValues */
        // Xprev[i] = X[i];
        // Yprev[i] = Y[i];

        /* .2 wurde erstzt durch calcNewValues*/
        // X[i] += VX[i] * tstep;
        // Y[i] += VY[i] * tstep;

        if((Xprev[i]>RoomXSize)&&(X[i]<=RoomXSize)) {
            NInRoom++;
        }
        if((Xprev[i]<=RoomXSize)&&(X[i]>RoomXSize)) {
            NInRoom--;

        }
    }


    /* .3 and .4 */
    for(i=0; i<N; i++) {

        /* (a) if the particle is on the board, its book-keeping
           arrays are modified only if its block has changed during
           the last update
           (b) if the particle is off-board, it will be removed */


        /* a */
        if(X[i]<XS) {
            j_old =   (int)floor(Xprev[i]*GX/XS)
                      + G*(int)floor(Yprev[i]*GY/YS);
            j_new = (int)floor(X[i]*GX/XS) + G*(int)floor(Y[i]*GY/YS);
            if( j_new != j_old ) {

                /* deleting particle i from its old block */
                j = j_old;
                if(BIndBd[j]==i) {
                    BIndBd[j] = BInd[i];
                } else {
                    j = BIndBd[j];
                    while(BInd[j]!=i) {
                        j = BInd[j];
                    }
                    BInd[j] = BInd[i];
                }


                /* inserting particle i into its new block */
                j = j_new;
                if(BIndBd[j]==-1) {
                    BIndBd[j] = i;
                    BInd[i] = -1;
                } else {
                    j = BIndBd[j];
                    while(BInd[j]!=-1) {
                        j = BInd[j];
                    }
                    BInd[j] = i;
                    BInd[i] = -1;
                }
            }
        } else {
            RemoveParticle( &N, i );
            i--;
        }
    }



    /* 4 */

    /* 4.1 */
    E[UpdNum+1] = 0.0;
    for(i=0; i<N; i++) {
        E[UpdNum+1] += VX[i] * cos(Phi[i]) + VY[i] * sin(Phi[i]);
    }
    if(N>0) {
        E[UpdNum+1] /= N;
    }


    /* 4.2 */
    for(i=0; i<N; i++) {
        VX[i] = vxnew[i];
        VY[i] = vynew[i];
        V[i] = sqrt(SQR(VX[i])+SQR(VY[i]));
        Vdir[i] = atan2(VY[i],VX[i]);
        Phi[i] = DirectionOfExit( i );
    }


    /* 4.3 */
    SimTime[UpdNum+1] = SimTime[UpdNum] + tstep;
    UpdNum++;


    /*
    if(NInjured>0){
    fprintf(stdout,"t[%d]=%g\n",UpdNum,SimTime[UpdNum]);
    fflush(stdout);
    }
    */


    /* 5 */
    free_vector(fwallx,0,allocN-1);
    free_vector(fwally,0,allocN-1);
    free_vector(fwpointx,0,allocN-1);
    free_vector(fwpointy,0,allocN-1);
    free_vector(fpairx,0,allocN-1);
    free_vector(fpairy,0,allocN-1);
    free_vector(fsmokex,0,allocN-1);
    free_vector(fsmokey,0,allocN-1);

    free_vector(fspx,0,allocN-1);
    free_vector(fspy,0,allocN-1);
    free_vector(fsumx,0,allocN-1);
    free_vector(fsumy,0,allocN-1);
    free_vector(vxnew,0,allocN-1);
    free_vector(vynew,0,allocN-1);
    free_vector(ftmagsum,0,allocN-1);
    free_vector(fcolx,0,allocN-1);
    free_vector(fcoly,0,allocN-1);


    /* 6 */
}
