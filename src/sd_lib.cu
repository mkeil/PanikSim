// library file for sd / panic used by sd.c and sd_crunch.c

int MainSwitch_DEFAULT = 0;
char *IFN_DEFAULT = "sd.par";
char *OFN_DEFAULT = "sd.dat";
char *OF2N_DEFAULT = "sd.dat2";

/********  global constants **********/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include "sd_lib.h"

#include "nrutil.c"
#include "base.c"
#include "readpar.c"
#include "Xlibext.c"



#define SD_LIB_EXIT {_E("sd_lib.c: Exiting to system.\n");exit(-1);}
#define MY_STRLEN 200
float EPSILON = 1.0e-5;

char *ParticleColorName[]= {"yellow", "grey"};
char *SmokeColorName="grey";


/* number of walls and wpoints */
#define NW 9
#define NWP 4
typedef struct wall {
    float x1, y1, x2, y2;
} wall;
typedef struct wpoint {
    float x, y;
} wpoint;


/********* global parameters -- to be read from parameter file **********/

static int N0, InjurySwitch, ColumnSwitch,
       X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim,
       SaveUN, DrawUN, Sleep, Draw,
       RndSeed, MaxUpdNum, AyS;

static float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean,
       deltaD, A, B, A_fire, B_fire, Kappa, C_Young, R, R_fire, V0, Tau,
       GaMe, GaTh, GaCM,
       SmokeStartTime, VSmoke, FCrush_over_1m,
       ColumnCenterX, ColumnCenterY, ColumnD,
       X11_Magn, SaveST, DrawST, DrawDMult, MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;
static char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN],
       X11_FontName[MY_STRLEN];


int *IPar[]= {&N0, &InjurySwitch, &ColumnSwitch,
              &X11_Margin, &X11_InFW,
              &X11_InFH, &X11_TLH, &X11_GrFH, &X11_RightRim, &SaveUN,
              &DrawUN, &Sleep, &Draw,
              &RndSeed, &MaxUpdNum, &AyS
             };
char *IParName[]= {"N0", "InjurySwitch",
                   "ColumnSwitch", "X11_Margin",
                   "X11_InFW", "X11_InFH", "X11_TLH", "X11_GrFH",
                   "X11_RightRim", "SaveUN", "DrawUN", "Sleep",
                   "Draw", "RndSeed", "MaxUpdNum", "AyS"
                  };
float *FPar[]= {&RoomXSize, &RoomYSize, &DoorWidth, &WallWidth, &Dmean,
                &deltaD, &A, &B, &A_fire, &B_fire, &Kappa,
                &C_Young, &R, &R_fire, &V0, &Tau,
                &GaMe, &GaTh, &GaCM, &SmokeStartTime, &VSmoke,
                &FCrush_over_1m,
                &ColumnCenterX, &ColumnCenterY, &ColumnD,
                &X11_Magn, &SaveST, &DrawST,
                &DrawDMult, &MaxSimTime, &Vmax, &H, &DefaultDeltaT,
                &C_NS, &V_ChangeLimit
               };
char *FParName[]= {"RoomXSize", "RoomYSize", "DoorWidth", "WallWidth",
                   "Dmean", "deltaD", "A", "B", "A_fire", "B_fire",
                   "Kappa", "C_Young",
                   "R", "R_fire", "V0", "Tau", "GaMe", "GaTh", "GaCM",
                   "SmokeStartTime", "VSmoke",
                   "FCrush_over_1m",
                   "ColumnCenterX", "ColumnCenterY", "ColumnD",
                   "X11_Magn", "SaveST",
                   "DrawST", "DrawDMult",
                   "MaxSimTime", "Vmax", "H",
                   "DefaultDeltaT", "C_NS", "V_ChangeLimit"
                  };
char *SPar[]= {BackGroundColorName, InfoColorName, X11_FontName};
char *SParName[]= {"BackGroundColorName", "InfoColorName",
                   "X11_FontName"
                  };


int IParNum = sizeof(IPar)/sizeof(int*),
    FParNum = sizeof(FPar)/sizeof(float*),
    SParNum = sizeof(SPar)/sizeof(char*);



/******* global variables ************/

int MainSwitch, UpdNum, N, GX, GY, G, Mb, Me, *BIndBd, *BInd, GaussFlag, NInRoom, *Injured, NInjured;
int X11_WWi, X11_WHe, BGCCode, ICCode, PaCNum, *PaCCode, SmokeCCode;

float *SimTime, XS, YS, *D, *Phi, *X, *Y, *Xprev, *Yprev, *V, *VX, *VY, *E, *Vdir, GaussSet1, GaussSet2, *V0of;
float FW_x;

XColor BGC_sdef,IC_sdef,PC_sdef[100],SmokeColor_sdef;

char IFN[MY_STRLEN], OFN[MY_STRLEN], OF2N[MY_STRLEN];
FILE *OFP,*OF2P;
wall *W;
wpoint *WP;

struct stat IFStatBuf;
long int IFModTime;


/********************* Funktionen *****************/

// void WallParticleRelation(int iw, int i, float *r, int *can_see)
// {
    // can_see: whether partice i is within the range of wall iw;  r: distance


    // switch(iw) {
    // case 0: {
            // *r = Y[i];
            // break;
        // }
    // case 1: {
            // *r = WP[0].x-X[i];
            // break;
        // }
    // case 2: {
            // *r = Y[i]-WP[0].y;
            // break;
        // }
    // case 3: {
            // *r = X[i]-WP[1].x;
            // break;
        // }
    // case 4: {
            // *r = RoomYSize-Y[i];
            // break;
        // }
    // case 5: {
            // *r = X[i]-WP[2].x;
            // break;
        // }
    // case 6: {
            // *r = WP[2].y-Y[i];
            // break;
        // }
    // case 7: {
            // *r = WP[3].x-X[i];
            // break;
        // }
    // case 8: {
            // *r = X[i];
            // break;
        // }
    // }


    // switch(iw) {
    // case 0: {
            // if(Y[i]<=R) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 1: {
            // if((X[i]>=WP[0].x-R)&&(X[i]<=WP[0].x)&&(Y[i]<=WP[0].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 2: {
            // if((X[i]>=WP[0].x)&&(X[i]<=WP[1].x)&&(Y[i]<=WP[0].y+R)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 3: {
            // if((X[i]>=WP[1].x)&&(X[i]<=WP[1].x+R)&&(Y[i]<=WP[1].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 4: {
            // if(Y[i]>=RoomYSize-R) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 5: {
            // if((X[i]>=WP[2].x)&&(X[i]<=WP[2].x+R)&&(Y[i]>=WP[2].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 6: {
            // if((X[i]>=WP[3].x)&&(X[i]<=WP[2].x)&&(Y[i]>=WP[2].y-R)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 7: {
            // if((X[i]<=WP[3].x)&&(X[i]>=WP[3].x-R)&&(Y[i]>=WP[3].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 8: {
            // if(X[i]<=R) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // }
// }
/*------------------------------*/

// void WallPsychForce(int iw, int i, float r, float *fx, float *fy)
// {

// #define tmp_f (A*exp(-(r-0.5*D[i])/B))

    // switch(iw) {
    // case 0: {
            // *fx = 0.0;
            // *fy = tmp_f;
            // break;
        // }
    // case 1: {
            // *fx = - tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 2: {
            // *fx = 0.0;
            // *fy = tmp_f;
            // break;
        // }
    // case 3: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 4: {
            // *fx = 0.0;
            // *fy = - tmp_f;
            // break;
        // }
    // case 5: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 6: {
            // *fx = 0.0;
            // *fy = - tmp_f;
            // break;
        // }
    // case 7: {
            // *fx = - tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 8: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // }

// #undef tmp_f
// }

// /*------------------------------*/

// void WallYoungForce(int iw, int i, float r, float *fx, float *fy)
// {

// #define tmp_f (2.0*C_Young*(0.5*D[i]-r))

    // switch(iw) {
    // case 0: {
            // *fx = 0.0;
            // *fy = tmp_f;
            // break;
        // }
    // case 1: {
            // *fx = - tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 2: {
            // *fx = 0.0;
            // *fy = tmp_f;
            // break;
        // }
    // case 3: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 4: {
            // *fx = 0.0;
            // *fy = - tmp_f;
            // break;
        // }
    // case 5: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 6: {
            // *fx = 0.0;
            // *fy = - tmp_f;
            // break;
        // }
    // case 7: {
            // *fx = - tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // case 8: {
            // *fx = tmp_f;
            // *fy = 0.0;
            // break;
        // }
    // }

// #undef tmp_f
// }


// /*------------------------------*/

// void WallTangForce_FS1( int iw, int i, float r, float *fx, float *fy )
// {

// #define tmp_delta_r (0.5*D[i]-r)

    // /* friction forces */
    // switch(iw) {
    // case 0:
    // case 2:
    // case 4:
    // case 6: {
            // *fx = -Kappa*tmp_delta_r*VX[i];
            // *fy = 0.0;
            // break;
        // }
    // case 1:
    // case 3:
    // case 5:
    // case 7:
    // case 8: {
            // *fx = 0.0;
            // *fy = -Kappa*tmp_delta_r*VY[i];
            // break;
        // }
    // }

// #undef tmp_delta_r
// }


/*------------------------------------*/

// void WPointParticleRelation(int iwp, int i, float *r, int *can_see)
// {
    // /* can_see: whether partice i is within the range of wpoint iwp
       // r: distance */

    // *r = sqrt(SQR(WP[iwp].x-X[i])+SQR(WP[iwp].y-Y[i]));

    // switch(iwp) {
    // case 0: {
            // if((X[i]<=WP[0].x)&&(Y[i]>=WP[0].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 1: {
            // if((X[i]>=WP[1].x)&&(Y[i]>=WP[1].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 2: {
            // if((X[i]>=WP[2].x)&&(Y[i]<=WP[2].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // case 3: {
            // if((X[i]<=WP[3].x)&&(Y[i]<=WP[3].y)) {
                // *can_see=1;
            // } else {
                // *can_see=0;
            // }
            // break;
        // }
    // }
// }

// /*------------------------------*/

// void WPointPsychForce(int iwp, int i, float r, float *fx, float *fy)
// {
    // /* exerted by wpoint iwp on particle i */

// #define tmp_f_over_r (A*exp(-(r-0.5*D[i])/B)/r)

    // *fx = (X[i]-WP[iwp].x) * tmp_f_over_r;
    // *fy = (Y[i]-WP[iwp].y) * tmp_f_over_r;

// #undef tmp_f_over_r
// }

// /*------------------------------*/

// void WPointYoungForce(int iwp, int i, float r, float *fx, float *fy)
// {
    // /* exerted by wpoint iwp on particle i */

    // float rx,ry;

// #define tmp_f_over_r ( 2.0*C_Young*(0.5*D[i]-r) / r)

    // rx=WP[iwp].x-X[i];
    // ry=WP[iwp].y-Y[i];
    // *fx = - rx * tmp_f_over_r;
    // *fy = - ry * tmp_f_over_r;

// #undef tmp_f_over_r
// }

// /*------------------------------*/

// void WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy)
// {
    // /* exerted by wpoint iwp on particle i */

    // float rx,ry,scal_prod_over_rsqr;

    // rx = X[i]-WP[iwp].x;
    // ry = Y[i]-WP[iwp].y;
    // scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(r);
    // *fx = -Kappa * (0.5*D[i]-r) * (   ry * scal_prod_over_rsqr );
    // *fy = -Kappa * (0.5*D[i]-r) * ( - rx * scal_prod_over_rsqr );
// }

/*------------------------------*/

void PP_PsychForce(int i1, int i2, float r, float *fx, float *fy)
{

    float f_over_r;

    f_over_r = A*exp(-(r-0.5*(D[i1]+D[i2]))/B) / r;
    *fx = (X[i1]-X[i2]) * f_over_r;
    *fy = (Y[i1]-Y[i2]) * f_over_r;
}

/*------------------------------*/

void PP_YoungForce(int i1, int i2, float r, float *fx, float *fy)
{

    float f_over_r;

    f_over_r = 2.0*C_Young*(0.5*(D[i1]+D[i2])-r) / r;
    *fx = (X[i1]-X[i2]) * f_over_r;
    *fy = (Y[i1]-Y[i2]) * f_over_r;
}

/*---------------------------*/

void PP_TangForce_FS1(int i1, int i2, float r, float *fx, float *fy)
{
    /* exerted by particle i2 on particle i1 */

    float rx,ry,vx,vy,scal_prod_over_rsqr;

    rx = X[i1]-X[i2];
    ry = Y[i1]-Y[i2];
    vx = VX[i1]-VX[i2];
    vy = VY[i1]-VY[i2];
    scal_prod_over_rsqr = (ry*vx - rx*vy) / SQR(r);
    *fx = -Kappa * (0.5*(D[i1]+D[i2])-r) * (   ry * scal_prod_over_rsqr );
    *fy = -Kappa * (0.5*(D[i1]+D[i2])-r) * ( - rx * scal_prod_over_rsqr );
}

/*---------------------------*/

float DirectionOfExit( int i )
{
    /* direction of exit for particle i */

    float dsqr, /* sqr of particle center - door-post distance */
          rsqr; /* sqr of particle's radius */


    /* behind the upper door-post */
    if((Y[i]<=0.5*YS-0.5*DoorWidth+0.5*D[i]+EPSILON)&&(X[i]<=RoomXSize)) {

        dsqr = SQR(W[1].x2-X[i]) + SQR(W[1].y2-Y[i]);
        rsqr = SQR(0.5*D[i])+EPSILON;
        if(dsqr<=rsqr) {
            /* very close to the door-post */
            if(Y[i]<=0.5*YS-0.5*DoorWidth) {
                return( 0.5*PI );
            } else {
                return(   0.5*PI
                          + atan2( W[1].y2-Y[i],W[1].x2-X[i] )
                      );
            }
        } else {
            /* well apart from the door-post */
            return(   atan2( 1.0, sqrt(dsqr/rsqr-1.0) )
                      + atan2( W[1].y2-Y[i],W[1].x2-X[i] )
                  );
        }
    }


    /* behind the lower door-post */
    else if((Y[i]>=0.5*YS+0.5*DoorWidth-0.5*D[i]-EPSILON)&&(X[i]<=RoomXSize)) {

        dsqr = SQR(W[6].x2-X[i]) + SQR(W[6].y2-Y[i]);
        rsqr = SQR(0.5*D[i])+EPSILON;
        if(dsqr<=rsqr) {
            /* very close to the door-post */
            if(Y[i]>=0.5*YS+0.5*DoorWidth) {
                return( -0.5*PI );
            } else {
                return( - 0.5*PI
                        + atan2( W[6].y2-Y[i],W[6].x2-X[i] )
                      );
            }
        } else {
            /* well apart from the door-post */
            return( - atan2( 1.0, sqrt(dsqr/rsqr-1.0) )
                    + atan2( W[6].y2-Y[i],W[6].x2-X[i] )
                  );
        }
    }


    /* in the center or outside */
    else {
        return 0.0;
    }
}

/*-----------------------------------------------*/

void RemoveParticle( int *n, int i )
{
    /* *n: number of particles now
       i: index of particle to be removed */

    int j;


    /* (a) particle i (which is off-board now)
       is removed from the book-keeping
       (block determined by previous coordinates)
       (b) if i != *n-1
           (b1) particle *n - 1 is removed from the book-keeping
          (block determined by previous coordinates)
           (b2) copying all values of particle *n-1 into i's place
           (b3) inserting particle i (that used to be indexed *n-1) into the
                book-keeping, into the block given by the previous
          coordinates (Xprev[i],Yprev[i]), and not into the block
          given by (X[i],Y[i])
          . reason: after this substitution (*n-1 -> i)
          particle i will be looked for in the block of
          (Xprev[i],Yprev[i]) in Upd
          because no one tells the main cycle (located in Upd),
          whether this particle is the result of a substitution
          or not
       (c) decrement particle number  ( *n to *n - 1 ) */


    /* a */
    j = (int)floor(Xprev[i]*GX/XS) + G*(int)floor(Yprev[i]*GY/YS);
    if(BIndBd[j]==i) {
        BIndBd[j] = BInd[i];
    } else {
        j = BIndBd[j];
        while(BInd[j]!=i) {
            j = BInd[j];
        }
        BInd[j] = BInd[i];
    }



    /* b */
    if(i!=*n-1) {

        /* b1 */

        j = (int)floor(Xprev[*n-1]*GX/XS) + G*(int)floor(Yprev[*n-1]*GY/YS);
        if(BIndBd[j]==*n-1) {
            BIndBd[j] = BInd[*n-1];
        } else {
            j = BIndBd[j];
            while(BInd[j]!=*n-1) {
                j = BInd[j];
            }
            BInd[j] = BInd[*n-1];
        }



        /* b2 */
        D[i] = D[*n-1];
        Phi[i] = Phi[*n-1];
        X[i] = X[*n-1];
        Y[i] = Y[*n-1];
        V[i] = V[*n-1];
        VX[i] = VX[*n-1];
        VY[i] = VY[*n-1];
        Xprev[i] = Xprev[*n-1];
        Yprev[i] = Yprev[*n-1];
        Vdir[i]=Vdir[*n-1];
        Injured[i]=Injured[*n-1];
        V0of[i]=V0of[*n-1];



        /* b3 */
        j = (int)floor(Xprev[i]*GX/XS) + G*(int)floor(Yprev[i]*GY/YS);
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




    /* c */
    (*n)--;
}

/*-----------------------------------------------------------------------*/

void EulTStep( float *tstep, float f )
{
    /* adjusts the time step in a way that the force (fx,fy) doesn't change the velocity of particle i by more than V_ChangeLimit */

    while( f*(*tstep) >= V_ChangeLimit ) {
        *tstep *= C_NS;
    }
}

/********************************/

void Start_Bare( int narg, char *argstr[] )
{
    /* reading command line parameters and parameter file */


    /* 1 */

    switch(narg) {
    case 1: {
            MainSwitch = MainSwitch_DEFAULT;
            strcpy(IFN,IFN_DEFAULT);
            strcpy(OFN,OFN_DEFAULT);
            strcpy(OF2N,OF2N_DEFAULT);
            break;
        }
    case 3: {
            MainSwitch = atoi(argstr[1]);
            strcpy(IFN,argstr[2]);
            strcpy(OFN,OFN_DEFAULT);
            strcpy(OF2N,OF2N_DEFAULT);
            break;
        }

    default: {
            _E("Usage:  either \"sd <MainSwitch> <input file name>\",");
            _E("        or     \"sd\" .\n");
            SD_LIB_EXIT;
        }
    }
    fprintf(stderr,"(default values are:\n MainSwitch = %d, input = %s, output = %s, %s)\n", MainSwitch,IFN,OFN,OF2N);
    fflush(stderr);

    /* 2 */

    /* reading parameters */
    readpar ( "start", IFN, IPar, IParName, IParNum,
              FPar, FParName, FParNum,
              SPar, SParName, SParNum );


    /* 3 */

}

/*----------------------------------------*/

void Init_Bare( char *init_switch )
{
    /* 1 global vars, mem.alloc.
       2 walls, wpoints
       3 particles
     */

    int i,j,ok_flag;



    /* 1 */
    stat(IFN, &IFStatBuf);
    IFModTime = IFStatBuf.st_mtime;
    UpdNum = 0;
    Mb = 0;
    Me = AyS-1;
    SimTime = vector( Mb, Me );
    SimTime[0] = 0.0;
    srand(RndSeed);

    XS = RoomXSize+WallWidth+X11_RightRim+EPSILON;
    YS = RoomYSize;
    N = N0;
    NInRoom = N0;

    GX = (int)MAX(1.0,floor(XS/R));
    GY = (int)MAX(1.0,floor(YS/R));
    G = (int)MAX(GX,GY);

    BIndBd = ivector( 0, SQR(G)-1 );
    BInd = ivector( 0, N0-1 );
    D = vector( 0, N0-1 );
    Phi = vector( 0, N0-1 );
    X = vector( 0, N0-1 );
    Y = vector( 0, N0-1 );
    Xprev = vector( 0, N0-1 );
    Yprev = vector( 0, N0-1 );
    V = vector( 0, N0-1 );
    VX = vector( 0, N0-1 );
    Vdir = vector(0,N0-1);
    VY = vector( 0, N0-1 );
    V0of = vector( 0, N0-1 );
    E = vector( Mb, Me );
    E[0] = 1.0;
    Injured = ivector(0,N0-1);



    /* 2 walls, wpoints */
    /* allocating memory:
     * if there's a column at the door,
     * the four faces and corners of the column have to be initialized,
     * too
     */
    W = (wall*)calloc(NW,sizeof(wall));
    WP = (wpoint*)calloc(NWP,sizeof(wpoint));



    /* every wall rotated by PI/2 points towards the inside of the room */

    /* upper part */
    W[0].x1 = 0.0;
    W[0].y1 = 0.0;
    W[0].x2 = XS;
    W[0].y2 = 0.0;

    W[1].x1 = RoomXSize;
    W[1].y1 = 0.0;
    W[1].x2 = W[2].x1 = WP[0].x = RoomXSize;
    W[1].y2 = W[2].y1 = WP[0].y = 0.5*RoomYSize-0.5*DoorWidth;
    W[2].x2 = W[3].x1 = WP[1].x = RoomXSize+WallWidth;
    W[2].y2 = W[3].y1 = WP[1].y = 0.5*RoomYSize-0.5*DoorWidth;
    W[3].x2 = RoomXSize+WallWidth;
    W[3].y2 = 0.0;


    /* lower part */
    W[4].x1 = XS;
    W[4].y1 = RoomYSize;
    W[4].x2 = 0.0;
    W[4].y2 = RoomYSize;

    W[5].x1 = RoomXSize+WallWidth;
    W[5].y1 = RoomYSize;

    W[5].x2 = W[6].x1 = WP[2].x = RoomXSize+WallWidth;
    W[5].y2 = W[6].y1 = WP[2].y = 0.5*RoomYSize+0.5*DoorWidth;
    W[6].x2 = W[7].x1 = WP[3].x = RoomXSize;
    W[6].y2 = W[7].y1 = WP[3].y = 0.5*RoomYSize+0.5*DoorWidth;
    W[7].x2 = RoomXSize;
    W[7].y2 = RoomYSize;


    /* left wall of the room */
    W[8].x1 = 0.0;
    W[8].y1 = RoomYSize;
    W[8].x2 = 0.0;
    W[8].y2 = 0.0;




    /* 3 */

    /* diameters and coordinates */
    for(i=0; i<N; i++) {
        D[i] =   (Dmean + deltaD)
                 - 2.0*deltaD * rand()/(RAND_MAX+1.0);
        X[i] =   0.5*H*D[i]+EPSILON
                 + (RoomXSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);
        Y[i] =   0.5*H*D[i]+EPSILON
                 + (RoomYSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);


        /* checking whether far enough from the column */
        ok_flag = 1;
        switch(ColumnSwitch) {
        default:
        case 0: {
                break;
            }
        case 1: {
                if(   SQR(X[i]-ColumnCenterX)+SQR(Y[i]-ColumnCenterY)
                        <= SQR(0.5*(D[i]+ColumnD))+EPSILON
                  ) {
                    ok_flag = 0;
                    i--;
                }
                break;
            }
        }


        /* checking distances to already existing particles */
        if(ok_flag==1) {
            for(j=0; j<i; j++) {
                if(     SQR(X[j] - X[i])
                        + SQR(Y[j] - Y[i])
                        <= SQR( 0.5*H*(D[i]+D[j]) ) + EPSILON
                  ) {
                    i = i - 1;
                    j = i - 1;
                }
            }
        }
    }


    /* book-keeping */
    for(i=0; i<SQR(G); i++) {
        BIndBd[i] = -1;
    }
    for(i=0; i<N; i++) {
        BInd[i] = -1;
    }

    for(i=0; i<N; i++) {
        j = (int)floor(X[i]*GX/XS) + G * (int)floor(Y[i]*GY/YS);

        if(BIndBd[j]==-1) {
            BIndBd[j] = i;
        } else {
            j = BIndBd[j];
            while(BInd[j]!=-1) {
                j = BInd[j];
            }
            BInd[j] = i;
        }
    }



    /* injuries, velocities and preferred directions */
    NInjured = 0;
    for(i=0; i<N; i++) {
        Injured[i] = 0;
    }

    for(i=0; i<N; i++) {
        Phi[i] = DirectionOfExit( i );
        Vdir[i] = Phi[i];
        V[i]=0.0;
        V0of[i]=V0;
        VX[i]=0.0;
        VY[i]=0.0;
    }
}

/*------------------------------*/

float GaussRand( float gmean, float gtheta, float gcutmult )
{
    /* generates a random number (x) with
       P(x) = exp[- (x-gmean)^2 / (2*gtheta)], if x is in
              [gmean - gcutmult*sqrt(gtheta), gmean + gcutmult*sqrt(gtheta)]
            = 0                              , if not */

    if( (GaussFlag==1) && (fabs(GaussSet2-gmean) <= gcutmult*sqrt(gtheta)) ) {
        GaussFlag = 0;
        return GaussSet2;
    } else {
        float v1,v2,rsq,fac;

        GaussFlag = 0;
        do {
            do {
                v1 = 1.0 - 2.0*(rand()/(RAND_MAX+1.0));
                v2 = 1.0 - 2.0*(rand()/(RAND_MAX+1.0));
            } while((rsq=v1*v1+v2*v2) >= 1.0);
            fac = sqrt(-2.0*gtheta*log(rsq)/rsq);
            GaussSet1 = v1*fac;
            GaussSet2 = v2*fac;
        } while(    (fabs(GaussSet1-gmean) > gcutmult*sqrt(gtheta))
                    && (fabs(GaussSet2-gmean) > gcutmult*sqrt(gtheta)) );

        if(fabs(GaussSet1-gmean) <= gcutmult*sqrt(gtheta)) {
            GaussFlag = 1;
            return GaussSet1;
        } else {
            GaussFlag = 0;
            return GaussSet2;
        }
    }
}

/*------------------------------*/

float EMean( char* sw, int unfreq, float stfreq )
{
    /* calculates the mean value of the efficiency of the system for the last
       few update steps -- NOTE: use this function only when UpdNum > 0

       if unfreq != 0, the average will be calculated for the last unfreq
       updates (the present one included)
       if unfreq == 0, the average will be calculated for the shortest
       possible time interval exceeding stfreq */

    int i, start;
    float e_mean, f;


    if(strcmp(sw,"un")==0) {
        start = UpdNum - unfreq;
    } else { /* i.e. if(strcmp(sw,"st")==0) */
        start = Mb; /* start from beginning of present time window */
        f = floor( SimTime[UpdNum] / stfreq );
        while( f - floor( SimTime[start] / stfreq ) > 1.0 ) {
            start++;
        }
        if( start==UpdNum ) {
            start--;
        }
    }
    e_mean = 0.0;
    for(i=start+1; i<=UpdNum; i++) {
        e_mean += E[i] * ( SimTime[i] - SimTime[i-1] );
    }
    e_mean /= SimTime[UpdNum] - SimTime[start];


    e_mean /= V0;
    return e_mean;
}


/*==============================*/



/*--------------------------------------------------*/

void Init_Demo()
{
    /* 1 general
       2 special
       */


    _E("Initializing, please wait... \n");


    /* 1 */
    Init_Bare("demo");



    /* 2 */

    /* opening files */
    if(!(OFP=fopen(OFN,"w"))) {
        fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",OFN);
        SD_LIB_EXIT;
    }
    fprintf(OFP,"UpdNum, SimTime, N, <E>\n");
    fflush(OFP);

    if(!(OF2P=fopen(OF2N,"w"))) {
        fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",OF2N);
        SD_LIB_EXIT;
    }
    fprintf(OF2P,"\n");
    fflush(OF2P);



    /* init visual or data output */
    X11_init();

    _E("... finished.\n");
}

/*------------------------------*/

void X11_init()
{
    int ii,last_ok;
    XColor sdef,edef;


    /* general */
    X11_WWi = X11_InFW + (int)(X11_Magn*XS) + X11_Margin;
    X11_WHe = (int)MAX( X11_InFH, X11_Magn*YS+X11_GrFH + 3*X11_Margin );
    g_win( "open", " self-driven", "sd", 0, 0, X11_WWi, X11_WHe, 4);
    g_font( "open", X11_FontName );


    /* colors */
    if( !XAllocNamedColor(display,cmap,BackGroundColorName,&edef,&sdef) ) {
        fprintf(stderr,"Error: couldn't allocate color: %s\n",
                BackGroundColorName);
        SD_LIB_EXIT;
    }
    BGCCode = sdef.pixel;

    if( !XAllocNamedColor(display,cmap,InfoColorName,&edef,&sdef) ) {
        fprintf(stderr,"Error: couldn't allocate color: %s\n",InfoColorName);
        SD_LIB_EXIT;
    }
    ICCode = sdef.pixel;


    PaCNum = sizeof(ParticleColorName)/sizeof(char*);
    PaCCode = ivector(0,PaCNum-1);
    for(ii=0,last_ok=0; ii<PaCNum; ii++) {
        if( !XAllocNamedColor(display,cmap,ParticleColorName[ii],&edef,&sdef) ) {
            fprintf(stderr,"WARNING: couldn't allocate color: %s\n",
                    ParticleColorName[ii]);
            fprintf(stderr,"Using %s instead\n",ParticleColorName[last_ok]);
            PaCCode[ii]=PaCCode[last_ok];
        } else {
            PaCCode[ii] = sdef.pixel;
            last_ok=ii;
        }
    }
}



/*------------------------------------------*/

void Pic()
{

    if(  (UpdNum==0)
            ||(  (UpdNum>0)
                 &&(  (  (DrawUN != 0)
                         &&(UpdNum % DrawUN == 0)
                      )
                      ||(  (DrawUN == 0)
                           &&(   floor( SimTime[UpdNum] / DrawST )
                                 > floor( SimTime[UpdNum-1] / DrawST )
                             )
                        )
                   )
              )
      ) {

        X11_Pic();

    }
}

/*------------------------------*/

void X11_Pic()
{
    /* 1 cleaning the whole window
       2 drawing particles (smoke front, column)
       3 walls
       4 cleaning the info surface, drawing info
       5 showing it, time delay
       */

    int i,disp_height;
    char disp_str[MY_STRLEN];
    /*  float pmean;*/
    float x;


    /* 1 */
    XSetForeground( display, gc, BGCCode );
    XFillRectangle( display, pix1, gc, 0, 0, X11_WWi, X11_WHe );


    /* 2 */
    for(i=0; i<N; i++) {
        XDrawParticle( i, X11_InFW, X11_Margin, X11_Magn);
    }



    /* 2.B */
    /* smoke front, if needed */
    if(  ((InjurySwitch==2)||(InjurySwitch==3))
            &&(SimTime[UpdNum]>=SmokeStartTime)) {
        XSetForeground( display, gc, ICCode );
        x =   X11_InFW + X11_Magn*(SimTime[UpdNum]-SmokeStartTime)*VSmoke;
        for(i=0; i<=X11_Magn*YS/6.0; i++) {
            XDrawLine(display, pix1, gc,
                      x, X11_Margin + 6*i,
                      x, X11_Margin + 6*i+3
                     );
        }
    }

    /* 2.C */
    /* column */
    switch(ColumnSwitch) {
    default:
    case 0: {
            break;
        }
    case 1: {
            XSetForeground( display, gc, ICCode );
            XDrawArc(display, pix1, gc,
                     (int)floor(X11_InFW+X11_Magn*(ColumnCenterX-0.5*ColumnD)),
                     (int)floor(X11_Margin+X11_Magn*(ColumnCenterY-0.5*ColumnD)),
                     (int)floor(X11_Magn*ColumnD),
                     (int)floor(X11_Magn*ColumnD),
                     0, 23040
                    );
            break;
        }
    }



    /* 3 */
    XSetForeground( display, gc, ICCode );
    for(i=0; i<NW; i++) {
        XDrawLine(display,pix1,gc,
                  (int)floor(X11_InFW+X11_Magn*W[i].x1),
                  (int)floor(X11_Margin+X11_Magn*W[i].y1),
                  (int)floor(X11_InFW+X11_Magn*W[i].x2),
                  (int)floor(X11_Margin+X11_Magn*W[i].y2)
                 );
    }





    /* 4 */
    XSetForeground( display, gc, BGCCode );

    /* cleaning the x=XS end of the field to allow particles
       leave the screen gradually */
    XFillRectangle( display, pix1, gc,
                    (int)floor(X11_InFW+X11_Magn*XS), 0,
                    (int)floor(X11_WWi-X11_InFW-X11_Magn*XS), X11_WHe );

    /* writing info */
    XSetForeground( display, gc, ICCode );
    disp_height = X11_Margin + X11_TLH;

    disp_height += X11_TLH;
    sprintf( disp_str, "t [%6d] = %.1f", UpdNum, SimTime[UpdNum] );
    XDrawString( display, pix1, gc, X11_Margin, disp_height,
                 disp_str, (signed int)strlen(disp_str) );

    disp_height += X11_TLH;
    sprintf( disp_str, "N (in room) = %d", NInRoom );
    XDrawString( display, pix1, gc, X11_Margin, disp_height,
                 disp_str, (signed int)strlen(disp_str) );

    disp_height += X11_TLH;
    sprintf( disp_str, "N_injured = %d", NInjured );
    XDrawString( display, pix1, gc, X11_Margin, disp_height,
                 disp_str, (signed int)strlen(disp_str) );

    disp_height += X11_TLH;
    sprintf( disp_str, "V0 = %g", V0 );
    XDrawString( display, pix1, gc, X11_Margin, disp_height,
                 disp_str, (signed int)strlen(disp_str) );

    disp_height += X11_TLH;
    sprintf( disp_str, "FWall_x = %.2f", FW_x );
    XDrawString( display, pix1, gc, X11_Margin, disp_height,
                 disp_str, (signed int)strlen(disp_str) );



    /* 5 */
    h_show(X11_WWi,X11_WHe);
    sleep(Sleep);
}


/*--------------------------------------*/



void XDrawParticle( int i, int leftxmargin, int upymargin, float magn)
{

    /* - drawing the particle  */

    int lxm = leftxmargin, uym = upymargin;
    float d,x,y;



    /* particle color */
    switch(Injured[i]) {
    case 0: {
            XSetForeground( display, gc, PaCCode[0] );
            break;
        }
    case 1: {
            XSetForeground( display, gc, PaCCode[1] );
            break;
        }
    }




    /* drawing particle */
    switch( Draw ) {
    default:
    case 0: {

            d = D[i];
            x = X[i];
            y = Y[i];

            XFillArc(display, pix1, gc,
                     (int)floor(lxm + magn * (x - d/2)),
                     (int)floor(uym + magn * (y - d/2)),
                     (int)floor(magn * d),
                     (int)floor(magn * d),
                     0, 23040
                    );

            break;
        }
    case 1: {

            d = D[i];
            x = X[i];
            y = Y[i];

            XFillArc(display, pix1, gc,
                     (int)floor(lxm + magn * (x - d/2)),
                     (int)floor(uym + magn * (y - d/2)),
                     (int)floor(magn * d),
                     (int)floor(magn * d),
                     0, 23040
                    );

            XSetForeground( display, gc, ICCode );
            XDrawArc(display, pix1, gc,
                     (int)floor(lxm + magn * (x - d/2)),
                     (int)floor(uym + magn * (y - d/2)),
                     (int)floor(magn * d),
                     (int)floor(magn * d),
                     0, 23040
                    );

            break;
        }
    case 2: {

            d = DrawDMult * D[i];
            x = X[i];
            y = Y[i];

            XFillArc(display, pix1, gc,
                     (int)floor(lxm + magn * (x - d/2)),
                     (int)floor(uym + magn * (y - d/2)),
                     (int)floor(magn * d),
                     (int)floor(magn * d),
                     0, 23040
                    );

            XSetForeground( display, gc, ICCode );
            XDrawArc(display, pix1, gc,
                     (int)floor(lxm + magn * (x - d/2)),
                     (int)floor(uym + magn * (y - d/2)),
                     (int)floor(magn * d),
                     (int)floor(magn * d),
                     0, 23040
                    );

            break;
        }
    case 3: {

            d = D[i];
            x = X[i];
            y = Y[i];

            if(sqrt(SQR(VX[i])+SQR(VY[i]))>=0.5) {
                XFillArc(display, pix1, gc,
                         (int)floor(lxm + magn * (x - d/2)),
                         (int)floor(uym + magn * (y - d/2)),
                         (int)floor(magn * d),
                         (int)floor(magn * d),
                         0, 23040
                        );
            }

            break;
        }
    }

}

/*-----------------------------------------------*/

void Save_Demo()
{

    char sw[MY_STRLEN];
    float simtime_now,simtime_now_minus_1,e_now,e_now_minus_1,
          e_mean;

    if(  (UpdNum==0)
            ||(  (UpdNum>0)
                 &&(  (  (SaveUN != 0)
                         &&(UpdNum % SaveUN == 0)
                      )
                      ||(  (SaveUN == 0)
                           &&(   floor( SimTime[UpdNum] / SaveST )
                                 > floor( SimTime[UpdNum-1] / SaveST )
                             )
                        )
                   )
              )
      ) {


        if(UpdNum>0) {
            if(SaveUN!=0) {
                strcpy(sw,"un");
            } else { /* ie. if(SaveUN==0) */
                strcpy(sw,"st");
            }
            e_mean=EMean(sw,SaveUN,SaveST);
        } else { /* ie. if(UpdNum==0) */
            e_mean=1.0;
        }
        fprintf(OFP,"%d\t%g\t%d\t%g\n",
                UpdNum,SimTime[UpdNum],N,e_mean);
        fflush(OFP);



        /* closing present time window, opening new time window */
        if(UpdNum>0) {

            simtime_now = SimTime[UpdNum];
            simtime_now_minus_1 = SimTime[UpdNum-1];
            e_now = E[UpdNum];
            e_now_minus_1 = E[UpdNum-1];

            free_vector(SimTime,Mb,Me);
            free_vector(E,Mb,Me);

            Mb = UpdNum-1;
            Me = UpdNum-1 + AyS-1;
            SimTime = vector(Mb,Me);
            SimTime[UpdNum-1] = simtime_now_minus_1;
            SimTime[UpdNum] = simtime_now;
            E = vector(Mb,Me);
            E[UpdNum-1] = e_now_minus_1;
            E[UpdNum] = e_now;
        }
    }
}

/*------------------------------*/

void Shutdown_Demo() {}

/*------------------------------*/

