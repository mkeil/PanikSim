#include "kernels.h"
#include "deviceFunc.cu"
#include "base.c"
#include <stdio.h>
__global__ void calcParticelForcesPar (float *fpairx, float *fpairy, float *ftmagsum, float *X, float *Y, float *D, float *VX, float *VY, int *Injured, int N, parameter *para)
{
    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;

    int m;
    float tmprsqr, tmpr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay;

    float R = para -> R;
    int InjurySwitch = para -> InjurySwitch;

    if (i < N) {
        for (m = 0; m < N; m ++) { 	// jedes Partikel überprüfen
            tmprsqr = SQR(X[i]-X[m]) + SQR(Y[i]-Y[m]);
            if( tmprsqr <= SQR(R) ) {
                tmpr = sqrt(tmprsqr);
                /* init */
                tmp_fpsx = tmp_fpsy = 0.0;
                tmp_fyox = tmp_fyoy = 0.0;
                tmp_ftax = tmp_ftay = 0.0;

                /* pair forces */
                /* Force(i,m,...) gives the force exerted by m on i, all forces are symmetric now */
                PP_PsychForce(i,m,tmpr,&tmp_fpsx,&tmp_fpsy, D, X, Y, para);
                if(tmpr<=0.5*(D[i]+D[m])) {
                    PP_YoungForce(i,m,tmpr,&tmp_fyox,&tmp_fyoy, D, X, Y, para);
                    PP_TangForce_FS1(i,m,tmpr,&tmp_ftax,&tmp_ftay, D, X, Y,VX, VY, para);
                }

                /* summing forces */
                if(Injured[i]==0) {
                    fpairx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
                    fpairy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
                } else { /* ie. if Injured[i]=1 */
                    fpairx[i] += tmp_fyox + tmp_ftax;
                    fpairy[i] += tmp_fyoy + tmp_ftay;
                }
                /* sum of magnitude of touching forces */
                if(InjurySwitch==1) {
                    ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
                }
            }
            if (i == m) { // verhindert dass die Werte verfälscht werden, wenn man mit sich selbst prüft
                ftmagsum[i] = 0.0;
                fpairx[i] = 0.0;
                fpairy[i] = 0.0;
            }

        }
    }
}






__global__ void calcParticelForces (float *fpairx, float *fpairy, float *ftmagsum, float *X, float *Y, float *D, float *VX, float *VY, int *Injured, float XS, float YS, int GX, int GY, int G, int *BIndBd, int *BInd, int N, parameter *para)
{
    // int b_ID = blockIdx.x;
    // int i =  b_ID * blockDim.x + threadIdx.x;

    int i;
    int j,k, l,mx, my, m;
    float tmprsqr, tmpr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay;

    float R = para -> R;
    int InjurySwitch = para -> InjurySwitch;
    for(i=0; i<N; i++) {
        // if (i < N) {

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
                                PP_PsychForce(i,m,tmpr,&tmp_fpsx,&tmp_fpsy, D, X, Y, para);
                                if(tmpr<=0.5*(D[i]+D[m])) {
                                    PP_YoungForce(i,m,tmpr,&tmp_fyox,&tmp_fyoy, D, X, Y, para);

                                    PP_TangForce_FS1(i,m,tmpr,&tmp_ftax,&tmp_ftay, D, X, Y,VX, VY, para);


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
}

__global__ void calcWallForces (float *fwallx, float *fwally, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nw)
{

    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;

    int iw;
    float R = para->R;
    int InjurySwitch = para->InjurySwitch;

    int can_see;
    float tmpr, tmp_fpsx, tmp_fpsy, tmp_fyox, tmp_fyoy, tmp_ftax, tmp_ftay ;

    if (i < N) {
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
    // if (i == 1) {
    // printf ("fwallx: &f", fwallx[1]);
    // }

}

__global__ void calcWPointForces (float *fwpointx, float *fwpointy, float *ftmagsum, float *D, int *Injured, float *X, float *Y, wpoint *WP, float *VX, float *VY, parameter *para, int N, int Nwp)
{

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

__global__ void calcColumnForces (float *fcolx, float *fcoly, float *ftmagsum,float *D, int *Injured, float *X, float *Y, float *VX, float *VY, parameter *para, int N)
{

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
    if (i < N) {

        switch(ColumnSwitch) {
        default:
        case 0: {

                fcolx[i] = fcoly[i] = 0.0;

                break;
            }
        case 1: {

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
                break;
            }
        }
    }
}

__global__ void calcInjuryForces (float *fsmokex, float *fsmokey, float *VX, float *VY, float *V0of, int *Injured, float *ftmagsum, int N, float SimTime, float *Phi, float *X, float *D, parameter *para)
{

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



    if (i < N) {


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

                if(SimTime>=SmokeStartTime) {
                    x_smokefront = (SimTime-SmokeStartTime)*VSmoke;


                    /* checking position compared to smoke front */
                    tmpr = X[i] - x_smokefront;

                    /* center of particle behind smoke front: injured */


                    if( tmpr < 0.5*D[i] ) {
                        if(Injured[i]==0) {
                            // printf("tmpr: %f, x_smokefront: %f \n", tmpr, x_smokefront);
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

__global__ void sumForces (float *fsumx,float *fsumy,  float *tStepVector, const float sqrt_fact, const float *VX,const float *VY,const float *V0of,const float *Phi,const float *fpairx,const float *fwallx,const float *fwpointx,const float *fpairy,const float *fwally,const float *fwpointy,const float *fsmokex,const float *fsmokey,const float *fcolx,const float *fcoly, const int N, const parameter *para)
{

    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;

    float Tau = para-> Tau;
    float DefaultDeltaT = para -> DefaultDeltaT;
    float V_ChangeLimit = para -> V_ChangeLimit;
    float C_NS = para -> C_NS;
    int InjurySwitch = para -> InjurySwitch;
    int ColumnSwitch = para -> ColumnSwitch;

    float fspx, fspy, ksi, eta;


    if (i < N) {

        /* self-propelling */
        fspx = 1/Tau * (V0of[i]*cos(Phi[i]) - VX[i]);
        fspy = 1/Tau * (V0of[i]*sin(Phi[i]) - VY[i]);



        // noise; die Verwendung habe ich erstmal rausgelassen, siehe erklärung in AnalyseSumForces
        // if(GaTh!=0.0) {
        // ksi = GaussRand(GaMe, GaTh, GaCM);
        // eta = 2.0*PI * rand() / (RAND_MAX+1.0);
        // } else {
        // ksi=0.0;
        // eta=0.0;
        // }

        ksi = 0.0;
        eta = 0.0;


        /* sum of forces */
        fsumx[i] =   fspx + fpairx[i] + fwallx[i] + fwpointx[i] + sqrt_fact * ksi * cos(eta);
        fsumy[i] =   fspy + fpairy[i] + fwally[i] + fwpointy[i] + sqrt_fact * ksi * sin(eta);

        // if (i == 1) {
        // printf ("Kraft in Summe Forces: %f, %f \n",fsumx[i], fsumy[i] );
        // }

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

        tStepVector[i] = EulTStep(DefaultDeltaT, sqrt(SQR(fsumx[i])+SQR(fsumy[i])), V_ChangeLimit, C_NS );

        // tStepVector[i] = 0.001;

        // tStepVector[i] = DefaultDeltaT;
        // float f = sqrt(SQR(fsumx[i])+SQR(fsumy[i]));
        // while ( f*(tStepVector[i]) >= V_ChangeLimit ) {
        // tStepVector[i] *= C_NS;

        // }


        // if ((i == 0) && (fsumx[0] < 0)) {
        //

        // }
        // printf ("fsumx: %f, fspx: %f, fpairx: %f, fwallx: %f, fwpointx: %f \n",fsumx[i], fspx, fpairx[i], fwallx[i], fwpointx[i] );
        // printf ("VX: %f, VY: %f \n", VX[0], VY[0]) ;
        // printf ("Phi[i]: %f \n ", Phi[i]);

    }
}


__global__ void NewVelocity (float *vxnew, float *vynew, const float *fsumx, const float *fsumy, const float *VX, const float *VY, const int *Injured, const int N, const float *tStepVector, parameter *para)
{

    /* new velocity */

    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;


    float vnew;
    float Vmax = para-> Vmax;
    int InjurySwitch = para -> InjurySwitch;

    if (i < N) {
        if(  (Injured[i]==1) &&((InjurySwitch==1)||(InjurySwitch==3))) {
            vxnew[i] = 0.0;
            vynew[i] = 0.0;
        } else {
            vxnew[i] = VX[i] + fsumx[i] * tStepVector[i];
            vynew[i] = VY[i] + fsumy[i] * tStepVector[i];
        }

        /* checking new velocity */
        vnew = sqrt( SQR(vxnew[i]) + SQR(vynew[i]) );
        if(vnew > Vmax) {
            vxnew[i] = vxnew[i]/vnew * Vmax;
            vynew[i] = vynew[i]/vnew * Vmax;
        }


    }
    // if (i == 0) {
    // printf ("Berechnung: alte Geschwindigkeit: %f, %f \n",VX[0], VY[0] );
    // printf ("Berechnung: Kraft: %f, %f \n",fsumx[0], fsumy[0] );
    // printf ("Berechnung: neue Geschwindigkeit: %f, %f \n",vxnew[0], vynew[0] );
    // }
}

__global__ void getNewValues (float* Xprev_d,float* X_d, float* Yprev_d, float *Y_d,float *VY_d, float* VX_d, int *NinRoomVektor_d, float tstep, int N, parameter *para)
{
    // speichert die alte Position und berechnet die neue Position der Partikel
    // es wird festgelegt, ob ein Partikel im Raum ist oder nicht
    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;

    float RoomXSize = para -> RoomXSize;

    if (i < N) {
        Xprev_d[i] = X_d[i];
        Yprev_d[i] = Y_d[i];
        X_d[i] += VX_d[i] * tstep;
        Y_d[i] += VY_d[i] * tstep;

        if((Xprev_d[i]>RoomXSize)&&(X_d[i]<=RoomXSize)) {
            NinRoomVektor_d[i] = 1;

        }
        if((Xprev_d[i] <=RoomXSize)&&(X_d[i]>RoomXSize)) {
            NinRoomVektor_d[i] = 0;

        }
        // if (i == 0) {
        // printf ("X_d : %f, Y_d: %f, VX: %f, VY: %f \n", X_d[i], Y_d[i], VX_d[i], VY_d[i]);
        // }

    }
}

__global__ void getMinTimeStep (const float *tStepVector, const int countElements, float *min)
{
    float erg = 10.0;
    int i;


    for (i = 0; i < countElements; i++) {

        if (tStepVector[i] < erg) {
            erg = tStepVector[i];
        }
    }

    *min = erg;
}


__global__ void sumUp (const int *summanden, const int countElements, int* sum)
{
    float summe = 0;
    int i;
    for (i = 0; i < countElements; i++) {
        summe = summe + summanden[i];
    }
    *sum = summe;
    // printf ("summe :%f \n", summe);
}

__global__ void storeNewVelocity (float *VX, float *VY, float *V, float *Vdir,  float *Phi, const float *X, const float *Y, const float *D, wall *W,  const float *vxnew, const float *vynew, parameter *para, const int N, float YS)
{

    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;
    // if (i == 0) {
    // printf ("Phi[i] old : %f \n", Phi[i]);
    // }

    if (i < N) {
        VX[i] = vxnew[i];
        VY[i] = vynew[i];
        V[i] = sqrt(SQR(VX[i])+SQR(VY[i]));
        Vdir[i] = atan2(VY[i],VX[i]);
        Phi[i] = DirectionOfExit(X[i], Y[i], D[i], YS, para, W);
    }
    // if (i == 0) {
    // printf ("Phi[i] new : %f \n", Phi[i]);
    // }
}

__global__ void setV0 (float *V0of, float V0, int N)
{
    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;

    if (i < N) {
        V0of[i] = V0;
    }
}

__global__ void setVdir_Phi (float *Vdir, float *Phi, int N, float *X, float *Y, float *D, float YS,parameter *para, wall *W)
{
    int b_ID = blockIdx.x;
    int i =  b_ID * blockDim.x + threadIdx.x;
    float dir;
    if (i < N) {
        dir = DirectionOfExit(X[i], Y[i], D[i], YS, para, W);
        Vdir[i] = dir;
        Phi[i] = dir;
    }
}

