#pragma once
#include "types.h"
#include "base.c"

void Save_Demo(int UpdNum, int  SaveUN, float SaveST, float *SimTime, float AyS, int Mb, int Me, float *E);

void calcParticelForcesOnHost (float *fpairx, float *fpairy, float *ftmagsum, float *X, float *Y, float *D, float *VX, float *VY, int *Injured, float XS, float YS, int GX, int GY, int G, int *BIndBd, int *BInd, int N, parameter *para);

void PP_PsychForceOnHost(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y, parameter *para);

void PP_YoungForceOnHost(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y, parameter *para);

void PP_TangForce_FS1OnHost(int i1, int i2, float r, float *fx, float *fy, float *D, float *X, float *Y,float *VX, float *VY, parameter *para);

void XDrawParticle(int leftxmargin, int upymargin, float magn, float d, float x, float y, int partikelInjured, X11props_t *X11props);

void X11_Pic(float XS, float YS, parameter *para, X11props_t *X11props, int N, int NInRoom, int NInjured, int UpdNum, float* SimTime, int NW, wall *W, float *D, float *X, float *Y, int *Injured);

bool needToDraw (int UpdNum, int  DrawUN, float DrawST, float *SimTime);

void X11_init(parameter *para, float XS, float YS, X11props_t *X11props);

void PlaceParticel (float *X, float *Y, float *D,parameter *para, float EPSILON);
	
void InitRoom (wall *W,wpoint *WP, parameter *para, float XS); 

void InitBookKeeping (int G, int N, float GX, float GY, float XS, float YS, int *BIndBd, int *BInd, float *X, float *Y); 

void prepareParameter (parameter* p);

void readpar ( char *sw, char *ifn, int *intValue[], char *intName[], int intNum, float *floatValue[], char *floatName[], int floatNum, char *stringValue[], char *stringName[], int stringNum );

void nrerror(char error_text[]);

float *vector(long nl, long nh);

int *ivector(long nl, long nh);

void free_vector(float *v, long nl, long nh); 

void free_ivector(int *v, long nl, long nh);