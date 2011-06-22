

#include "update.h"
#include <glog/logging.h>
#include <stdio.h>

__global__ void storeOldValues (float* Xprev_d,float* X_d, float* Yprev_d, float*Y_d)
{

    int i = threadIdx.x;
    // printf ("%d\n", i);
    if (i <= 200) {
        Xprev_d[i] = X_d[i];
        Yprev_d[i] = Y_d[i];
    }
}

__global__ void calcNewValues (float* X_d,float* Y_d,float* VY_d,float* VX_d, float tstep)
{


    int i = threadIdx.x;
    if (i <= 200) {
        X_d[i] += VX_d[i] * tstep;
        Y_d[i] += VY_d[i] * tstep;
    }

}
void Upd()
{

    int allocN,i,j,k,l,mx,my,m,can_see,iwp,iw,j_old,j_new;
    float *fwallx,*fwally,*fwpointx,*fwpointy,*fpairx,*fpairy,
          *fspx,*fspy,*fsumx,*fsumy,*vxnew,*vynew,tstep,tmpr,
          tmp_fpsx,tmp_fpsy,tmp_fyox,tmp_fyoy,tmp_ftax,tmp_ftay,
          tmprsqr,sqrt_fact,ksi,eta,vnew,*ftmagsum,*fsmokex,*fsmokey,
          x_smokefront,tmpf,f_over_r,scal_prod_over_rsqr,rx,ry,
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






    /* 1.1 */
    /* wall force */
    for(i=0; i<N; i++) {
        for(iw=0; iw<NW; iw++) {

            WallParticleRelation(iw,i,&tmpr,&can_see);
            if((can_see==1)&&(tmpr<=R)) {

                /* init */
                tmp_fpsx = tmp_fpsy = 0.0;
                tmp_fyox = tmp_fyoy = 0.0;
                tmp_ftax = tmp_ftay = 0.0;

                /* psychological force */
                WallPsychForce(iw,i,tmpr,&tmp_fpsx,&tmp_fpsy);
                /* Young and tangential forces */
                if(tmpr<=0.5*D[i]) {
                    WallYoungForce(iw,i,tmpr,&tmp_fyox,&tmp_fyoy);

                    WallTangForce_FS1(iw,i,tmpr,&tmp_ftax,&tmp_ftay);

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

                /* measuring x component of touching force exerted
                on walls left and right from exit
                 -- only in demo mode */

                if((iw==1)||(iw==7)) {
                    FW_x -= tmp_fyox + tmp_ftax;
                }
            }
        }
    }



    /* 1.2 */
    /* wpoint force */
    for(i=0; i<N; i++) {
        for(iwp=0; iwp<NWP; iwp++) {

            WPointParticleRelation(iwp,i,&tmpr,&can_see);
            if((can_see==1)&&(tmpr<=R)) {

                /* init */
                tmp_fpsx = tmp_fpsy = 0.0;
                tmp_fyox = tmp_fyoy = 0.0;
                tmp_ftax = tmp_ftay = 0.0;

                /* computing forces */
                WPointPsychForce(iwp,i,tmpr,&tmp_fpsx,&tmp_fpsy);
                if(tmpr<=0.5*D[i]) {
                    WPointYoungForce(iwp,i,tmpr,&tmp_fyox,&tmp_fyoy);

                    WPointTangForce_FS1(iwp,i,tmpr,&tmp_ftax,&tmp_ftay);

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

                /* measuring x component of touching force exerted
                   on walls left and right from exit
                   -- only in demo mode */

                if((iwp==0)||(iwp==3)) {
                    FW_x -= tmp_fyox + tmp_ftax;
                }

            }
        }
    }



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



    /* 1.4
     * column
     */
    switch(ColumnSwitch) {
    default:
    case 0: {
            for(i=0; i<N; i++) {
                fcolx[i] = fcoly[i] = 0.0;
            }
            break;
        }
    case 1: {
            for(i=0; i<N; i++) {
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
            }
            break;
        }
    }



    /* 1.5 */
    /* injuries */

    switch(InjurySwitch) {
    case 0: {
            break;
        }
    case 1: {

            /* case: people crushed */
            for(i=0; i<N; i++) {

                /* newly injured */
                if((ftmagsum[i]>FCrush_over_1m*PI*D[i])&&(Injured[i]==0)) {
                    Injured[i] = 1;
                    NInjured++;
                    V0of[i] = 0.0;
                }
            }
            break;
        }
    case 2:
    case 3: {

            /* case: smoke front */
            if(SimTime[UpdNum]>=SmokeStartTime) {
                x_smokefront = (SimTime[UpdNum]-SmokeStartTime)*VSmoke;

                for(i=0; i<N; i++) {
                    /* checking position compared to smoke front */
                    tmpr = X[i] - x_smokefront;

                    /* center of particle behind smoke front: injured */
                    if( tmpr < 0.5*D[i] ) {
                        if(Injured[i]==0) {
                            Injured[i] = 1;
                            NInjured++;
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
            }
            break;
        }
    }



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
    float *Xprev_d, *Yprev_d, *X_d, *Y_d, *VX_d, *VY_d;
    int sizeFloatVector = N * sizeof(float);

    cudaError_t error = cudaMalloc (&Xprev_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


    error = cudaMalloc (&Yprev_d, sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


    error = cudaMalloc (&X_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


    error = cudaMalloc (&Y_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&VY_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&VX_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    // alte Werte hochkopieren
    error = cudaMemcpy(X_d, X, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(Y_d, Y, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(VY_d, VY, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMemcpy(VX_d, VX, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


    storeOldValues <<<1, 200>>> (Xprev_d,X_d,Yprev_d,Y_d);
    cudaThreadSynchronize();


    calcNewValues <<<1, N>>> (X_d,Y_d,VY_d,VX_d,tstep);
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
