#include "hostFunc.h"

#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#include <glog/logging.h>


#define NR_END 1
#define FREE_ARG char*
#define READPAR_EXIT {fprintf(stderr,"readpar EXIT\n");fflush(stderr);exit(-1);}
#define SD_LIB_EXIT {_E("sd_lib.c: Exiting to system.\n");exit(-1);}

#include "Xlib_mod.h"

#define MY_STRLEN 200


void XDrawParticle(int leftxmargin, int upymargin, float magn, float d, float x, float y, int partikelInjured, X11props_t *X11props)
{
    // nötige Props zur Verfügung stellen

    Display *display = X11props -> display;
    GC gc = X11props -> gc;
    Pixmap pix1 = X11props -> pix1;
    int *PaCCode = X11props -> PaCCode;




    /* - drawing the particle  */
    int lxm = leftxmargin, uym = upymargin;

    /* particle color */
    switch(partikelInjured) {
    case 0: {
            XSetForeground( display, gc, PaCCode[0] );
            break;
        }
    case 1: {
            XSetForeground( display, gc, PaCCode[1] );
            break;
        }
    }

    XFillArc(display, pix1, gc,
             (int)floor(lxm + magn * (x - d/2)),
             (int)floor(uym + magn * (y - d/2)),
             (int)floor(magn * d),
             (int)floor(magn * d),
             0, 23040
            );
}
void X11_Pic(float XS, float YS, parameter *para, X11props_t *X11props, int N, int NInRoom, int NInjured, int UpdNum, float* SimTime, int NW, wall *W, float *D, float *X, float *Y, int *Injured)

{

    // Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    GC gc = X11props -> gc;
    Pixmap pix1 = X11props -> pix1;
    int BGCCode = X11props -> BGCCode;
    int ICCode = X11props -> ICCode;
    int X11_WWi = X11props -> X11_WWi;
    int X11_WHe = X11props -> X11_WHe;

    // parameter auslesen
    int X11_Margin = para -> X11_Margin;
    float X11_Magn = para -> X11_Magn;
    int X11_TLH = para -> X11_TLH;
    int Sleep = para -> Sleep;
    int InjurySwitch = para -> InjurySwitch;
    int ColumnSwitch = para -> ColumnSwitch;
    float V0 = para -> V0;
    float SmokeStartTime = para -> SmokeStartTime;
    float VSmoke = para -> VSmoke;
    int X11_InFW = para -> X11_InFW;
    float ColumnCenterX = para -> ColumnCenterX;
    float ColumnD = para -> ColumnD;
    float ColumnCenterY = para -> ColumnCenterY;


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
        XDrawParticle(X11_InFW, X11_Margin, X11_Magn, D[i], X[i], Y[i], Injured[i], X11props);
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

    // cleaning the x=XS end of the field to allow particles leave the screen gradually
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
    sprintf( disp_str, "N = %d", N );
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




    /* 5 */
    h_show(X11_WWi,X11_WHe, X11props);
    sleep(Sleep);
}


bool needToDraw (int UpdNum, int  DrawUN, float DrawST, float *SimTime)
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

        return true;

    }
    return false;
}



void Save_Demo(int UpdNum, int  SaveUN, float SaveST, float *SimTime,  float AyS, int Mb, int Me, float *E)
{
    LOG (INFO) << "Save Demo gestartet. ";
    float simtime_now,simtime_now_minus_1,e_now,e_now_minus_1;

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




        // closing present time window, opening new time window
        if(UpdNum>0) {
            LOG (INFO) << "closing present time window, opening new time window";
            simtime_now = SimTime[UpdNum];
            simtime_now_minus_1 = SimTime[UpdNum-1];
            e_now = E[UpdNum];
            e_now_minus_1 = E[UpdNum-1];
            LOG (INFO) << "Werte berechnet.";

            free_vector(SimTime,Mb,Me);
            free_vector(E,Mb,Me);
            LOG (INFO) << "Vektoren freigegeben";

            Mb = UpdNum-1;
            Me = UpdNum-1 + AyS-1;
            SimTime = vector(Mb,Me);
            SimTime[UpdNum-1] = simtime_now_minus_1;
            SimTime[UpdNum] = simtime_now;
            E = vector(Mb,Me);
            E[UpdNum-1] = e_now_minus_1;
            E[UpdNum] = e_now;
            LOG (INFO) << "neue Vektoren erzeugt.";
        }
    }
    LOG (INFO) << "Save Demo beendet.";
}

void X11_init(parameter *para, float XS, float YS, X11props_t *X11props)
{
    LOG (INFO) << "x11 init gestartet.";
    int ii,last_ok;
    XColor sdef,edef;

    // benötige Proberties lokal zur Verfügung stellen
    Display *display = X11props -> display;
    Colormap cmap = X11props -> cmap;



    LOG (INFO) << "Props lokal zur Verfügung gestellt.";

    // parameter auslesen
    int X11_InFW = para -> X11_InFW;
    int X11_InFH = para -> X11_InFH;
    int X11_Margin = para -> X11_Margin;
    int X11_GrFH = para -> X11_GrFH;

    float X11_Magn = para -> X11_Magn;

    char X11_FontName [MY_STRLEN];
    strcpy (X11_FontName, para->X11_FontName);
    char BackGroundColorName [MY_STRLEN];
    strcpy (BackGroundColorName, para -> BackGroundColorName);
    char InfoColorName[MY_STRLEN];
    strcpy (InfoColorName, para -> InfoColorName);


    LOG (INFO) << "Parameter ausgelesen.";

    char *ParticleColorName[]= {"yellow", "grey"};
    

    int BGCCode, ICCode, PaCNum, *PaCCode;
    LOG (INFO) << "lokale Werte initialisiert." ;

    /* general */
    int X11_WWi = X11_InFW + (int)(X11_Magn*XS) + X11_Margin;
    int X11_WHe = (int)MAX( X11_InFH, X11_Magn*YS+X11_GrFH + 3*X11_Margin );

    g_win( "open", " PanicSimulator", "PanSim", 0, 0, X11_WWi, X11_WHe, 4, X11props);
    LOG (INFO) << "g_win durchgeführt.";
    g_font( "open", X11_FontName, X11props);
    LOG (INFO) << "g_font durchgeführt.";

    // nach den Aufrufen von g_win und g_font müßen display und cmap neu gesetzt werden
    display = X11props -> display;
    cmap = X11props -> cmap;

    /* colors */
    if( !XAllocNamedColor(display,cmap,BackGroundColorName,&edef,&sdef) ) {
        fprintf(stderr,"Error: couldn't allocate color: %s\n", BackGroundColorName);
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

    // lokale Veränderungen speichern
    X11props -> display = display;
    X11props -> cmap = cmap ;

	
	
    X11props -> BGCCode = BGCCode;
    X11props -> ICCode = ICCode;
    X11props -> PaCNum = PaCNum;
    X11props -> PaCCode = PaCCode;
    
    X11props -> X11_WWi = X11_WWi;
    X11props -> X11_WHe = X11_WHe;
}

void PlaceParticel (float *X, float *Y, float *D,parameter *para, float EPSILON)
{

    int i,j,ok_flag;
    int N = para-> N0;
    int ColumnSwitch = para -> ColumnSwitch;
    float ColumnCenterX = para -> ColumnCenterX;
    float ColumnCenterY = para -> ColumnCenterY;
    float ColumnD = para -> ColumnD;
    float Dmean = para-> Dmean;
    float deltaD = para -> deltaD;
    float H = para -> H;
    float RoomXSize = para -> RoomXSize;
    float RoomYSize = para -> RoomYSize;


    for(i=0; i<N; i++) {
        D[i] =   (Dmean + deltaD) - 2.0*deltaD * rand()/(RAND_MAX+1.0);
        X[i] =   0.5*H*D[i]+EPSILON + (RoomXSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);
        Y[i] =   0.5*H*D[i]+EPSILON + (RoomYSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);


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
}


void InitRoom (wall *W,wpoint *WP, parameter *para, float XS)
{

    float RoomXSize = para -> RoomXSize;
    float RoomYSize = para -> RoomYSize;
    float DoorWidth = para -> DoorWidth;
    float WallWidth = para -> WallWidth;



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
}

void InitBookKeeping (int G, int N, float GX, float GY, float XS, float YS, int *BIndBd, int *BInd, float *X, float *Y)
{
    int i, j;
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
}






void prepareParameter (parameter* p)
{
    LOG(INFO) << "pepareParamter aufgerufen";
    // die Parameter werden eingelesen und in die passende Struct eingefügt
    int N0, InjurySwitch, ColumnSwitch, X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim, SaveUN, DrawUN, Sleep, Draw, RndSeed, MaxUpdNum, AyS;

    float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean, deltaD, A, B, A_fire, B_fire, Kappa, C_Young, R, R_fire, V0, Tau, GaMe, GaTh, GaCM, SmokeStartTime, VSmoke, FCrush_over_1m, ColumnCenterX, ColumnCenterY, ColumnD, X11_Magn, SaveST, DrawST, DrawDMult, MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;

    char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN], X11_FontName[MY_STRLEN];
    LOG(INFO) << "Variablen deklariert.";

    int *IPar[]= {&N0, &InjurySwitch, &ColumnSwitch, &X11_Margin, &X11_InFW, &X11_InFH, &X11_TLH, &X11_GrFH, &X11_RightRim, &SaveUN, &DrawUN, &Sleep, &Draw, &RndSeed, &MaxUpdNum, &AyS };
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
    char *SParName[]= {"BackGroundColorName", "InfoColorName", "X11_FontName"};

    LOG(INFO) << "Arrays intialisiert.";

    int IParNum = sizeof(IPar)/sizeof(int*),
        FParNum = sizeof(FPar)/sizeof(float*),
        SParNum = sizeof(SPar)/sizeof(char*);

    LOG(INFO) << "ReadPar wird aufgerufen.";

    readpar ( "start", "panic.par", IPar, IParName, IParNum, FPar, FParName, FParNum, SPar, SParName, SParNum );

    LOG(INFO) << "ReadPar fertig.";

    // paras in Struct kopieren
    // int Para;
    p->N0	=	N0	;
    p->InjurySwitch	=	 InjurySwitch	;
    p->ColumnSwitch	=	 ColumnSwitch	;
    p->X11_Margin	=	 X11_Margin	;
    p->X11_InFW	=	 X11_InFW	;
    p->X11_InFH	=	 X11_InFH	;
    p->X11_TLH	=	 X11_TLH	;
    p->X11_GrFH	=	 X11_GrFH	;
    p->X11_RightRim	=	 X11_RightRim	;
    p->SaveUN	=	 SaveUN	;
    p->DrawUN	=	 DrawUN	;
    p->Sleep	=	 Sleep	;
    p->Draw	=	 Draw	;
    p->RndSeed	=	 RndSeed	;
    p->MaxUpdNum	=	 MaxUpdNum	;
    p->AyS	=	 AyS	;
    LOG(INFO) << "IntPar kopiert.";
    // float Parameter
    p->RoomXSize	=	RoomXSize	;
    p->RoomYSize	=	 RoomYSize	;
    p->DoorWidth	=	 DoorWidth	;
    p->WallWidth	=	 WallWidth	;
    p->Dmean	=	 Dmean	;
    p->deltaD	=	 deltaD	;
    p->A	=	 A	;
    p->B	=	 B	;
    p->A_fire	=	 A_fire	;
    p->B_fire	=	 B_fire	;
    p->Kappa	=	 Kappa	;
    p->C_Young	=	 C_Young	;
    p->R	=	 R	;
    p->R_fire	=	 R_fire	;
    p->V0	=	 V0	;
    p->Tau	=	 Tau	;
    p->GaMe	=	 GaMe	;
    p->GaTh	=	 GaTh	;
    p->GaCM	=	 GaCM	;
    p->SmokeStartTime	=	 SmokeStartTime	;
    p->VSmoke	=	 VSmoke	;
    p->FCrush_over_1m	=	 FCrush_over_1m	;
    p->ColumnCenterX	=	 ColumnCenterX	;
    p->ColumnCenterY	=	 ColumnCenterY	;
    p->ColumnD	=	 ColumnD	;
    p->X11_Magn	=	 X11_Magn	;
    p->SaveST	=	 SaveST	;
    p->DrawST	=	 DrawST	;
    p->DrawDMult	=	 DrawDMult	;
    p->MaxSimTime	=	 MaxSimTime	;
    p->Vmax	=	 Vmax	;
    p->H	=	 H	;
    p->DefaultDeltaT	=	 DefaultDeltaT	;
    p->C_NS	=	 C_NS	;
    p->V_ChangeLimit	=	 V_ChangeLimit	;
    LOG(INFO) << "Float Parameter kopiert.";

    // char parameter
    strcpy (p->BackGroundColorName, BackGroundColorName);
    strcpy (p->InfoColorName, InfoColorName);
    strcpy (p->X11_FontName, X11_FontName);
    LOG(INFO) << "Char Parameter kopiert.";

}

/* readpar_v_2000_03_02.c

   reading parameter file with this format:
   1st column     #
   2nd            one of these characters: i f s (integer/float/string)
   3rd            one of these characters: 0 1 (1:interactive parameter,0:not)
   4th            parameter name
   5th            parameter value
*/


void readpar ( char *sw, char *ifn,
               int *intValue[], char *intName[], int intNum,
               float *floatValue[], char *floatName[], int floatNum,
               char *stringValue[], char *stringName[], int stringNum )
{

    /* sw: switch = "start" or "re"
       ifn: input file name
       each array starts with the 0. element
       */


    FILE *ifp;
    int ii,i,*intFound,*floatFound,*stringFound,exitFlag,
        interactive,tmpInt;
    char tmpString[100], tmpStringVal[100],c;
    float tmpFloat;


    /* malloc, init, etc */
    intFound = ivector(0,intNum-1);
    floatFound = ivector(0,floatNum-1);
    stringFound = ivector(0,stringNum-1);
    for(i=0; i<intNum; i++) {
        intFound[i]=0;
    }
    for(i=0; i<floatNum; i++) {
        floatFound[i]=0;
    }
    for(i=0; i<stringNum; i++) {
        stringFound[i]=0;
    }

    if( !(ifp = fopen(ifn, "r")) ) {
        fprintf(stderr,"readpar: Couldn't open \"%s\" for reading\n",ifn);
        READPAR_EXIT;
    }


    for(ii=0; ii<intNum+floatNum+stringNum; ii++) {
        while(getc(ifp)!=0x23) {}; /* 0x23 = '#' */

        /* type of parameter */
        fscanf(ifp,"%s",&c);
        switch(c) {
        case 0x69: { /* 0x69 = 'i' */
                exitFlag=0;
                fscanf(ifp,"%d %s %d",&interactive,tmpString,&tmpInt);
                for(i=0; (exitFlag==0)&&(i<intNum); i++) {
                    if(strcmp(intName[i],tmpString)==0) {
                        intFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            *(intValue[i])=tmpInt;
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this integer parameter: %s\n",tmpString);
                }
                break;
            }
        case 0x66: { /* 0x66 = 'f' */
                exitFlag=0;
                fscanf(ifp,"%d %s %f",&interactive,tmpString,&tmpFloat);
                for(i=0; (exitFlag==0)&&(i<floatNum); i++) {
                    if(strcmp(floatName[i],tmpString)==0) {
                        floatFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            *(floatValue[i])=tmpFloat;
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this float parameter: %s\n",tmpString);
                }
                break;
            }
        case 0x73: { /* 0x73 = 's' */
                exitFlag=0;
                fscanf(ifp,"%d %s %s",&interactive,tmpString,tmpStringVal);
                for(i=0; (exitFlag==0)&&(i<stringNum); i++) {
                    if(strcmp(stringName[i],tmpString)==0) {
                        stringFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            strcpy(stringValue[i],tmpStringVal);
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this string parameter: %s\n",tmpString);
                }
                break;
            }
        }
    }


    /* checking whether all parameters have been found */
    for(i=0; i<intNum; i++) {
        if(intFound[i]==0) {
            fprintf(stderr,"readpar ERROR: integer parameter %s not found in %s\n",intName[i],ifn);
            READPAR_EXIT;
        }
    }
    for(i=0; i<floatNum; i++) {
        if(floatFound[i]==0) {
            fprintf(stderr,"readpar ERROR: float parameter %s not found in %s\n",floatName[i],ifn);
            READPAR_EXIT;
        }
    }
    for(i=0; i<stringNum; i++) {
        if(stringFound[i]==0) {
            fprintf(stderr,"readpar ERROR: stringparameter %s not found in %s\n",stringName[i],ifn);
            READPAR_EXIT;
        }
    }

    free_ivector(intFound,0,intNum-1);
    free_ivector(floatFound,0,floatNum-1);
    free_ivector(stringFound,0,stringNum-1);
}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
    fprintf(stderr,"Numerical Recipes run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
    float *v;

    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    if (!v) {
        nrerror("allocation failure in vector()");
    }
    return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
    int *v;

    v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
    if (!v) {
        nrerror("allocation failure in ivector()");
    }
    return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}