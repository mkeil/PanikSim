
#include "prepareParameter.h"
#define MY_STRLEN 200
#include "readpar.c"
void prepareParameter (parameter* p) {
	// die Parameter werden eingelesen und in die passende Struct eingefügt
	static int N0, InjurySwitch, ColumnSwitch, X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim, SaveUN, DrawUN, Sleep, Draw, RndSeed, MaxUpdNum, AyS;

	static float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean, deltaD, A, B, A_fire, B_fire, Kappa, C_Young, R, R_fire, V0, Tau, GaMe, GaTh, GaCM, SmokeStartTime, VSmoke, FCrush_over_1m, ColumnCenterX, ColumnCenterY, ColumnD, X11_Magn, SaveST, DrawST, DrawDMult, MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;

	static char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN], X11_FontName[MY_STRLEN];

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

	int IParNum = sizeof(IPar)/sizeof(int*),
		FParNum = sizeof(FPar)/sizeof(float*),
		SParNum = sizeof(SPar)/sizeof(char*);
	
	readpar ( "start", "sd.par", IPar, IParName, IParNum, FPar, FParName, FParNum, SPar, SParName, SParNum );
	
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
	
	// char parameter
	strcpy (p->BackGroundColorName, BackGroundColorName);
	strcpy (p->InfoColorName, InfoColorName);
	strcpy (p->X11_FontName, X11_FontName);
}