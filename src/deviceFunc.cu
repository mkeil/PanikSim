#include "deviceFunc.h"
#include "base.c"

#include <stdio.h>

__device__ void PP_TangForce_FS1(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y,float *VX, float *VY, parameter *para)
{
    /* exerted by particle i2 on particle i1 */
    float Kappa = para -> Kappa;
    float rx,ry,vx,vy,scal_prod_over_rsqr;

    rx = X[i1]-X[i2];
    ry = Y[i1]-Y[i2];
    vx = VX[i1]-VX[i2];
    vy = VY[i1]-VY[i2];
    scal_prod_over_rsqr = (ry*vx - rx*vy) / SQR(r);
    *fx = -Kappa * (0.5*(D[i1]+D[i2])-r) * (   ry * scal_prod_over_rsqr );
    *fy = -Kappa * (0.5*(D[i1]+D[i2])-r) * ( - rx * scal_prod_over_rsqr );
}

__device__ void PP_PsychForce(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y, parameter *para)
{

    float A = para -> A;
    float B = para -> B;

    float f_over_r;

    f_over_r = A*exp(-(r-0.5*(D[i1]+D[i2]))/B) / r;
    *fx = (X[i1]-X[i2]) * f_over_r;
    *fy = (Y[i1]-Y[i2]) * f_over_r;
}

__device__ void PP_YoungForce(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y, parameter *para)
{
    float C_Young = para -> C_Young;
    float f_over_r;

    f_over_r = 2.0*C_Young*(0.5*(D[i1]+D[i2])-r) / r;
    *fx = (X[i1]-X[i2]) * f_over_r;
    *fy = (Y[i1]-Y[i2]) * f_over_r;
}


__device__ void WallParticleRelation(int iw, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP, parameter *para)
{
    // can_see: whether partice i is within the range of wall iw;  r: distance

    float RoomYSize = para-> RoomYSize;
    float R = para -> R;

    switch(iw) {
    case 0: {
            *r = yCoor;
            break;
        }
    case 1: {
            *r = WP[0].x -xCoor;
            break;
        }
    case 2: {
            *r = yCoor-WP[0].y;
            break;
        }
    case 3: {
            *r = xCoor-WP[1].x;
            break;
        }
    case 4: {
            *r = RoomYSize-yCoor;
            break;
        }
    case 5: {
            *r = xCoor-WP[2].x;
            break;
        }
    case 6: {
            *r = WP[2].y-yCoor;
            break;
        }
    case 7: {
            *r = WP[3].x-xCoor;
            break;
        }
    case 8: {
            *r = xCoor;
            break;
        }
    }


    switch(iw) {
    case 0: {
            if(yCoor<=R) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 1: {
            if((xCoor>=WP[0].x-R)&&(xCoor<=WP[0].x)&&(yCoor<=WP[0].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 2: {
            if((xCoor>=WP[0].x)&&(xCoor<=WP[1].x)&&(yCoor<=WP[0].y+R)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 3: {
            if((xCoor>=WP[1].x)&&(xCoor<=WP[1].x+R)&&(yCoor<=WP[1].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 4: {
            if(yCoor>=RoomYSize-R) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 5: {
            if((xCoor>=WP[2].x)&&(xCoor<=WP[2].x+R)&&(yCoor>=WP[2].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 6: {
            if((xCoor>=WP[3].x)&&(xCoor<=WP[2].x)&&(yCoor>=WP[2].y-R)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 7: {
            if((xCoor<=WP[3].x)&&(xCoor>=WP[3].x-R)&&(yCoor>=WP[3].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 8: {
            if(xCoor<=R) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    }
}




__device__ void WallTangForce_FS1( int iw, int i, float r, float *fx, float *fy, float diameter, float VelocityX_Dir, float VelocityY_Dir, parameter *para )
{
    float Kappa = para-> Kappa;
#define tmp_delta_r (0.5*diameter-r)

    /* friction forces */
    switch(iw) {
    case 0:
    case 2:
    case 4:
    case 6: {
            *fx = -Kappa*tmp_delta_r*VelocityX_Dir;
            *fy = 0.0;
            break;
        }
    case 1:
    case 3:
    case 5:
    case 7:
    case 8: {
            *fx = 0.0;
            *fy = -Kappa*tmp_delta_r*VelocityY_Dir;
            break;
        }
    }

#undef tmp_delta_r
}

__device__ void WallPsychForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para)
{
    float A = para-> A;
    float B = para -> B;
#define tmp_f (A*exp(-(r-0.5*diameter)/B))

    switch(iw) {
    case 0: {
            *fx = 0.0;
            *fy = tmp_f;
            break;
        }
    case 1: {
            *fx = - tmp_f;
            *fy = 0.0;
            break;
        }
    case 2: {
            *fx = 0.0;
            *fy = tmp_f;
            break;
        }
    case 3: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    case 4: {
            *fx = 0.0;
            *fy = - tmp_f;
            break;
        }
    case 5: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    case 6: {
            *fx = 0.0;
            *fy = - tmp_f;
            break;
        }
    case 7: {
            *fx = - tmp_f;
            *fy = 0.0;
            break;
        }
    case 8: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    }

#undef tmp_f
}

__device__ void WallYoungForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para)
{

    float C_Young = para-> C_Young;
#define tmp_f (2.0*C_Young*(0.5*diameter-r))

    switch(iw) {
    case 0: {
            *fx = 0.0;
            *fy = tmp_f;
            break;
        }
    case 1: {
            *fx = - tmp_f;
            *fy = 0.0;
            break;
        }
    case 2: {
            *fx = 0.0;
            *fy = tmp_f;
            break;
        }
    case 3: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    case 4: {
            *fx = 0.0;
            *fy = - tmp_f;
            break;
        }
    case 5: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    case 6: {
            *fx = 0.0;
            *fy = - tmp_f;
            break;
        }
    case 7: {
            *fx = - tmp_f;
            *fy = 0.0;
            break;
        }
    case 8: {
            *fx = tmp_f;
            *fy = 0.0;
            break;
        }
    }

#undef tmp_f
}

__device__ void WPointYoungForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para )

{
    /* exerted by wpoint iwp on particle i */
    float C_Young = para-> C_Young;
    float rx,ry;

#define tmp_f_over_r ( 2.0*C_Young*(0.5*diameter-r) / r)

    rx=WP.x-xCoor;
    ry=WP.y-yCoor;
    *fx = - rx * tmp_f_over_r;
    *fy = - ry * tmp_f_over_r;

#undef tmp_f_over_r
}

__device__ void WPointPsychForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para )

{
    /* exerted by wpoint iwp on particle i */
    float A = para-> A;
    float B = para-> B;
#define tmp_f_over_r (A*exp(-(r-0.5*diameter)/B)/r)

    *fx = (xCoor-WP.x) * tmp_f_over_r;
    *fy = (yCoor-WP.y) * tmp_f_over_r;

#undef tmp_f_over_r
}


__device__ void WPointParticleRelation(int iwp, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP)
{
    /* can_see: whether partice i is within the range of wpoint iwp r: distance */

    *r = sqrt(SQR(WP[iwp].x-xCoor)+SQR(WP[iwp].y-yCoor));

    switch(iwp) {
    case 0: {
            if((xCoor<=WP[0].x)&&(yCoor>=WP[0].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 1: {
            if((xCoor>=WP[1].x)&&(yCoor>=WP[1].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 2: {
            if((xCoor>=WP[2].x)&&(yCoor<=WP[2].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    case 3: {
            if((xCoor<=WP[3].x)&&(yCoor<=WP[3].y)) {
                *can_see=1;
            } else {
                *can_see=0;
            }
            break;
        }
    }
}

__device__ void WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, float VelocityX_Dir, float VelocityY_Dir, wpoint WP, parameter *para )

{
    /* exerted by wpoint iwp on particle i */
    float Kappa = para-> Kappa;
    float rx,ry,scal_prod_over_rsqr;

    rx = xCoor-WP.x;
    ry = yCoor-WP.y;
    scal_prod_over_rsqr = (ry*VelocityX_Dir - rx*VelocityY_Dir) / SQR(r);
    *fx = -Kappa * (0.5*diameter-r) * (   ry * scal_prod_over_rsqr );
    *fy = -Kappa * (0.5*diameter-r) * ( - rx * scal_prod_over_rsqr );
}

__device__ float EulTStep(float tmpTimeStep, float f, float V_ChangeLimit, float C_NS )
{
    /* adjusts the time step in a way that the force (fx,fy) doesn't change the velocity of particle i by more than V_ChangeLimit */

    while( f*(tmpTimeStep) >= V_ChangeLimit ) {
        tmpTimeStep *= C_NS;
    }

    return tmpTimeStep;
}


__device__ float DirectionOfExit(float xCoor, float yCoor, float diameter, float YS, parameter *para, wall *W)
{
    float DoorWidth = para-> DoorWidth;
    float RoomXSize = para -> RoomXSize;
    float EPSILON = 1.0e-5;

    // printf ("Dir of Exit: x: %f, y: %f, d: %f, YS: %f \n", xCoor, yCoor, diameter, YS);
    /* direction of exit for particle i */

    float dsqr, /* sqr of particle center - door-post distance */
          rsqr; /* sqr of particle's radius */


    /* behind the upper door-post */
    if((yCoor<=0.5*YS-0.5*DoorWidth+0.5*diameter+EPSILON)&&(xCoor<=RoomXSize)) {

        dsqr = SQR(W[1].x2-xCoor) + SQR(W[1].y2-yCoor);
        rsqr = SQR(0.5*diameter)+EPSILON;
        if(dsqr<=rsqr) {
            /* very close to the door-post */
            if(yCoor<=0.5*YS-0.5*DoorWidth) {
                return( 0.5*PI );
            } else {
                return(   0.5*PI
                          + atan2( W[1].y2-yCoor,W[1].x2-xCoor )
                      );
            }
        } else {
            /* well apart from the door-post */
            return(   atan2( 1.0, sqrt(dsqr/rsqr-1.0) )
                      + atan2( W[1].y2-yCoor,W[1].x2-xCoor )
                  );
        }
    }


    /* behind the lower door-post */
    else if((yCoor>=0.5*YS+0.5*DoorWidth-0.5*diameter-EPSILON)&&(xCoor<=RoomXSize)) {

        dsqr = SQR(W[6].x2-xCoor) + SQR(W[6].y2-yCoor);
        rsqr = SQR(0.5*diameter)+EPSILON;
        if(dsqr<=rsqr) {
            /* very close to the door-post */
            if(yCoor>=0.5*YS+0.5*DoorWidth) {
                return( -0.5*PI );
            } else {
                return( - 0.5*PI
                        + atan2( W[6].y2-yCoor,W[6].x2-xCoor )
                      );
            }
        } else {
            /* well apart from the door-post */
            return( - atan2( 1.0, sqrt(dsqr/rsqr-1.0) )
                    + atan2( W[6].y2-yCoor,W[6].x2-xCoor )
                  );
        }
    }


    /* in the center or outside */
    else {
        return 0.0;
    }
}