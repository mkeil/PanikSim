#pragma once

void Upd(parameter para_h);

void __device__ WallParticleRelation(int iw, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP, parameter *para);
void __device__ WallTangForce_FS1( int iw, int i, float r, float *fx, float *fy, float diameter, float VelocityX_Dir, float VelocityY_Dir, parameter *para );
void __device__ WallPsychForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para);
void __device__ WallYoungForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para);

void __device__ WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, float VelocityX_Dir, float VelocityY_Dir, wpoint WP, parameter *para );
void __device__ WPointYoungForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para );
void __device__ WPointPsychForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para );
void __device__ WPointParticleRelation(int iw, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP);