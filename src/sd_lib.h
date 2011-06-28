#pragma once


void PP_PsychForce(int i1, int i2, float r, float *fx, float *fy);
void PP_YoungForce(int i1, int i2, float r, float *fx, float *fy);
void PP_TangForce_FS1(int i1, int i2, float r, float *fx, float *fy);

float DirectionOfExit( int i );
void RemoveParticle( int *n, int i );
void EulTStep( float *tstep, float f );
void Start_Bare( int narg, char *argstr[] );
void Init_Bare( char *init_switch );
void Upd();
float GaussRand( float gmean, float gtheta, float gcutmult );
float EMean( char* sw, int cnfreq, float stfreq );
void Start_Demo( int narg, char *argstr[] );
void Init_Demo();
void X11_init();
void Pic();
void X11_Pic();
void XDrawParticle( int i, int leftxmargin, int upymargin, float magn);
void Save_Demo();
void Shutdown_Demo();
void ReInit();