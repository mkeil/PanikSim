
__device__ void WallParticleRelation(int iw, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP, parameter *para); 

__device__ void WallTangForce_FS1( int iw, int i, float r, float *fx, float *fy, float diameter, float VelocityX_Dir, float VelocityY_Dir, parameter *para ); 

__device__ void WallPsychForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para); 

__device__ void WallYoungForce(int iw, int i, float r, float *fx, float *fy, float diameter, parameter *para); 

__device__ void WPointYoungForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para ); 

__device__ void WPointPsychForce(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, wpoint WP, parameter *para );

__device__ void WPointParticleRelation(int iwp, int i, float *r, int *can_see, float yCoor, float xCoor, wpoint *WP); 

__device__ void WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy, float xCoor, float yCoor, float diameter, float VelocityX_Dir, float VelocityY_Dir, wpoint WP, parameter *para ); 

__device__ float EulTStep(float tmpTimeStep, float f, float V_ChangeLimit, float C_NS );

__device__ float DirectionOfExit(float xCoor, float yCoor, float diameter, float YS, parameter *para, wall *W); 