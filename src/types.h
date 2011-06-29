#pragma once
#define MY_STRLEN 200
struct parameter
{
    float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean, deltaD, A, B, A_fire, B_fire, Kappa, C_Young, R, R_fire, V0, Tau, GaMe, GaTh, GaCM, SmokeStartTime, VSmoke, FCrush_over_1m, ColumnCenterX, ColumnCenterY, ColumnD, X11_Magn, SaveST, DrawST, DrawDMult, MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;
	int N0, InjurySwitch, ColumnSwitch, X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim, SaveUN, DrawUN, Sleep, Draw, RndSeed, MaxUpdNum, AyS;
	char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN], X11_FontName[MY_STRLEN];
};

typedef struct parameter parameter;

/* number of walls and wpoints */
#define NW 9
#define NWP 4
typedef struct wall {
    float x1, y1, x2, y2;
} wall;
typedef struct wpoint {
    float x, y;
} wpoint;