#pragma once

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>

#define MY_STRLEN 200

struct parameter
{
    float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean, deltaD, A, B, A_fire, B_fire, Kappa, C_Young, R, R_fire, V0, Tau, GaMe, GaTh, GaCM, SmokeStartTime, VSmoke, FCrush_over_1m, ColumnCenterX, ColumnCenterY, ColumnD, X11_Magn, SaveST, DrawST, DrawDMult, MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;
	int N0, InjurySwitch, ColumnSwitch, X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim, SaveUN, DrawUN, Sleep, Draw, RndSeed, MaxUpdNum, AyS;
	char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN], X11_FontName[MY_STRLEN];
};

typedef struct parameter parameter;

struct X11props_t
{
	Display *display;
	Window win;
	GC gc;
	Pixmap pix1;              /* hidden windows */
	int xx,yy,wwidth,hheight;
	unsigned int bborder_width;
	Colormap cmap;
	Visual *vis;
	int screen;
	int BGCCode, ICCode, PaCNum, *PaCCode, SmokeCCode, X11_WWi, X11_WHe; // globale Variablen für die Darstellung der Simulation
};

typedef struct X11props_t X11props_t;


typedef struct wall {
    float x1, y1, x2, y2;
} wall;
typedef struct wpoint {
    float x, y;
} wpoint;