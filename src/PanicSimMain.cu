#include <glog/logging.h>
#include <stdio.h>

#include "types.h"
#include "hostFunc.h"
#include "kernels.h"

#include "base.c"




int main( int argc, char **argv)
{
    google::InitGoogleLogging(argv[0]);

    // A) Parameter einlesen
    parameter *para_h;
    para_h = (parameter*) malloc (sizeof(parameter));
    prepareParameter (para_h);
    LOG (INFO) << "Parameter eingelsen.";

    // globale Werte setzen
    int N0 = para_h-> N0;
    float EPSILON = 1.0e-5;
    float AyS = para_h -> AyS;
    int UpdNum = 0;
    int Mb = 0;
    int Me = AyS-1;
    float *SimTime = vector( Mb, Me );
    SimTime[0] = 0.0;
    srand(para_h->RndSeed);
    float XS = (para_h->RoomXSize)+(para_h->WallWidth)+(para_h->X11_RightRim)+EPSILON;
    float YS = para_h->RoomYSize;
    int N = N0;
    int NInRoom = N0;
    int NInjured = 0;
    int GX = (int)MAX(1.0,floor(XS/para_h->R));
    int GY = (int)MAX(1.0,floor(YS/para_h->R));
    int G = (int)MAX(GX,GY);
    float MaxSimTime = para_h -> MaxSimTime;
    int MaxUpdNum = para_h -> MaxUpdNum;
    float DefaultDeltaT = para_h -> DefaultDeltaT;

    int DrawUN = para_h -> DrawUN;
    float DrawST = para_h -> DrawST;

    float sqrt_fact;

    const int NW = 9; // Anzahl der Wände
    const int NWP = 4; // Anzahl der wpoints (Ecken)

    int *BIndBd = ivector( 0, SQR(G)-1 );
    int *BInd = ivector( 0, N0-1 );
    float *D = vector( 0, N0-1 );
    float *Phi = vector( 0, N0-1 );
    float *X = vector( 0, N0-1 );
    float *Y = vector( 0, N0-1 );
    float *Xprev = vector( 0, N0-1 );
    float *Yprev = vector( 0, N0-1 );
    float *V = vector( 0, N0-1 );
    float *VX = vector( 0, N0-1 );
    float *Vdir = vector(0,N0-1);
    float *VY = vector( 0, N0-1 );
    float *V0of = vector( 0, N0-1 );
    float *E = vector( Mb, Me );
    E[0] = 1.0;
    int *Injured = ivector(0,N0-1);
    LOG (INFO) << "globale Werte gesetzt.";

    dim3 dimBlock(256);
    dim3 dimGrid((N0 + dimBlock.x - 1) / dimBlock.x);

    LOG(INFO) << "Block dimensions: " << dimBlock.x << " " << dimBlock.y << " " << dimBlock.z;
    LOG(INFO) << "Grid dimensions: " << dimGrid.x << " " << dimGrid.y << " " << dimGrid.z;


    /* 2 walls, wpoints initalisieren */
    wall *W;
    wpoint *WP;
    W = (wall*)calloc(NW,sizeof(wall));
    WP = (wpoint*)calloc(NWP,sizeof(wpoint));
    InitRoom (W,WP, para_h, XS);
    LOG (INFO) << "Wände initialisiert.";

    // partikel plazieren
    PlaceParticel (X, Y, D,para_h, EPSILON);
    LOG (INFO) << "Partikel erzeugt & plaziert.";

    InitBookKeeping (G, N0, GX, GY, XS, YS, BIndBd, BInd, X, Y);
    LOG (INFO) << "Initiales Bookkeeping durchgeführt.";




    // alle device Pointer erzeugen;
    float *Xprev_d, *Yprev_d, *X_d, *Y_d, *V_d, *VX_d, *VY_d, *Vdir_d, *fwallx_d, *fwally_d,*fwpointx_d, *fwpointy_d, *ftmagsum_d, *D_d, *fcolx_d, *fcoly_d, *fsmokex_d, *fsmokey_d, *V0of_d, *SimTime_d, *Phi_d, *fsumx_d, *fsumy_d, *tStepVector_d, *fpairx_d, *fpairy_d, *timeStep_d, *vxnew_d, *vynew_d;

    int *BIndBd_d, *BInd_d;
    int sizeBIndBd_d = SQR (G) -1;



    int *Injured_d, *NInjured_d;

    int *NInRoom_d, *NinRoomVektor_d;

    wpoint *WP_d;

    wall *W_d;

    parameter *para_d;

    int sizeFloatVector = N0 * sizeof(float);

    cudaError_t error;

    // DeviceVektoren der Partikeleigenschaften
    error = cudaMalloc (&BInd_d,N0 * sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&BIndBd_d,sizeBIndBd_d * sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&Injured_d,N0 * sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&D_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&vxnew_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&vynew_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&V0of_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&Xprev_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&Yprev_d, sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&X_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&Y_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&VY_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&VX_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&V_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&Vdir_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&Phi_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    // DevicePointer der globalen Werte
    error = cudaMalloc (&timeStep_d,sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // TimeStep auf dem Device
    error = cudaMalloc (&NInRoom_d,sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // Anzahl der Personen im Raum auf dem Device
    error = cudaMalloc (&NInjured_d,sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // Anzahl der Verletzten
    error = cudaMalloc (&tStepVector_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // Hilfsvektor um den neuen minimalen Timestep zu finden
    error = cudaMalloc (&NinRoomVektor_d,N0 * sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error); // Hilfsvektor um den die neue Anzahl der Personen im Raum zu ermitteln

    error = cudaMalloc (&SimTime_d, sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    // DeviceVektoren für die Szene
    error = cudaMalloc (&W_d,NW * sizeof(wall));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&WP_d,NWP * sizeof(wpoint));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&para_d,sizeof(parameter));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    // Device Vektoren der Teilkräfte
    error = cudaMalloc (&fwallx_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fwally_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&fwpointx_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fwpointy_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&fpairx_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fpairy_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&fcolx_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fcoly_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&fsmokex_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fsmokey_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&fsumy_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMalloc (&fsumx_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    error = cudaMalloc (&ftmagsum_d,sizeFloatVector);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

    LOG (INFO) << "Device Pointer alloziert.";


    // nötige Elemente hochkopieren u.A. die Parameter auf Device kopieren  *****
    // Vektoren für die Szene hochkopieren

    error = cudaMemcpy(W_d, W, NW * sizeof(wall), cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(WP_d, WP, NWP * sizeof(wpoint), cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(para_d, para_h, sizeof(parameter), cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "Parameter und Wände hochkopiert.";

    // Eigenschaftsvektoren kopieren bzw. auf Device initialisieren

    // die folgenden leer initialiseren
    error = cudaMemset (Injured_d, 0, N * sizeof(int));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (vxnew_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (vynew_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (Xprev_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (Yprev_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (VY_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (VX_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (V_d, 0, N * sizeof(float));
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "Leere Eigenschaftsvektoren wurden initialisiert.";

    // erzeugte Werte hochkopieren
    error = cudaMemcpy(X_d, X, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(Y_d, Y, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(D_d, D, sizeFloatVector, cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(BInd_d, BInd, N * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemcpy(BIndBd_d, BIndBd, sizeBIndBd_d * sizeof(int), cudaMemcpyHostToDevice);
    CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);



    LOG (INFO) << "x Koordinaten, y Koordinaten und Durchmesser wurden kopiert.";


    // StartWert für V0 setzen
    float V0 = para_h-> V0;
    setV0 <<<dimGrid, dimBlock>>> (V0of_d, V0, N0);
    cudaThreadSynchronize(); error = cudaGetLastError (); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "Das Ergebniss von setV0 ist: " << cudaGetErrorString(error);


    // StartWerte für Vdir und Phi setzten --> beide mit direction of Exit initialisieren
    setVdir_Phi <<<dimGrid, dimBlock>>> (Vdir_d, Phi_d, N, X_d, Y_d, D_d, YS,para_d, W_d);
    cudaThreadSynchronize(); error = cudaGetLastError (); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "Das Ergebniss von setVdir_Phi ist: " << cudaGetErrorString(error);

    // DevicePointer der globalen Werte mit initialen Werten setzen
    error = cudaMemset (NInRoom_d, 0, sizeof(int)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (NInjured_d, 0, sizeof(int)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (timeStep_d, 0, sizeof(float)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (tStepVector_d, 0, N * sizeof(float)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (NinRoomVektor_d, 0, N * sizeof(int)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset (SimTime_d, 0, N * sizeof(float)); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "DevicePointer der globalen Werte mit initialen Werten gesetzt.";

    // temporäre Kraftverktoren initialisiert
    error = cudaMemset(fpairy_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fpairx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fwallx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fwally_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fcolx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fcoly_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fsmokex_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fsmokey_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fwpointx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fwpointy_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(ftmagsum_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fsumx_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    error = cudaMemset(fsumy_d, 0, sizeFloatVector); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
    LOG (INFO) << "temporäre Kraftverktoren initialisiert.";


    // X11 Props intialisieren und X11 Window init
    X11props_t *X11props;
    int sizeX11Porps = (int) sizeof(X11props_t);
    X11props = (X11props_t*) malloc (sizeX11Porps);
    // enthaltene Pointer initalsieren
    X11props->display = NULL;
    X11props->vis = NULL;

    LOG (INFO) << "sizeX11Porps: " << sizeX11Porps;
    X11_init(para_h, XS, YS, X11props);
    LOG (INFO) << "X11 Fenster initialisiert";



    
    float tstep;
    LOG (WARNING) << "Schleife gestartet.";
    do {
        
        LOG (INFO) << "UpDate Nr: " << UpdNum;
        error = cudaMemset(fpairy_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fpairx_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fwallx_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fwally_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fcolx_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fcoly_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fsmokex_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fsmokey_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fwpointx_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fwpointy_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(ftmagsum_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fsumx_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemset(fsumy_d, 0, sizeFloatVector);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        LOG (INFO) << "temporäre Vektoren initialisiert.";


        // Werte zum zeichnen runterkopieren
        error = cudaMemcpy(X, X_d, sizeFloatVector, cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemcpy(Y, Y_d, sizeFloatVector, cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemcpy(D, D_d, sizeFloatVector, cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemcpy(Injured, Injured_d, N0 * sizeof(int), cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemcpy(VX, VX_d, sizeFloatVector, cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        error = cudaMemcpy(VY, VY_d, sizeFloatVector, cudaMemcpyDeviceToHost);
        CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);


        LOG (INFO) << "Werte zum Zeichnen zurückkopiert.";

 


        if (needToDraw (UpdNum, DrawUN, DrawST, SimTime)) {
            X11_Pic(XS, YS, para_h, X11props, N, NInRoom, NInjured, UpdNum, SimTime, NW, W, D, X, Y, Injured);
            // printf (", gezeichnet");
        }
        LOG (INFO) << "Fertig mit zeichnen." ;


        // Update Schritt
        tstep = DefaultDeltaT;
        calcWallForces<<<dimGrid, dimBlock>>> (fwallx_d, fwally_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, WP_d, VX_d, VY_d, para_d, N, NW); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        calcWPointForces <<<dimGrid, dimBlock>>> (fwpointx_d, fwpointy_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, WP_d, VX_d, VY_d, para_d, N, NWP); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        LOG (INFO) << "W und WP Forces berechnet. Start CalcPartikelForces" ;

        //  particle-particle forces berechnen
        
        calcParticelForcesPar <<<dimGrid, dimBlock>>> (fpairx_d, fpairy_d, ftmagsum_d,X_d, Y_d, D_d, VX_d, VY_d, Injured_d, N, para_d);        
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        LOG (INFO) << "calcParticelForcesOnHost fertig" ;


        calcColumnForces <<<dimGrid, dimBlock>>> (fcolx_d, fcoly_d, ftmagsum_d, D_d, Injured_d, X_d, Y_d, VX_d, VY_d, para_d, N); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        calcInjuryForces <<<dimGrid, dimBlock>>> (fsmokex_d, fsmokey_d, VX_d, VY_d, V0of_d, Injured_d,ftmagsum_d, N, SimTime[UpdNum], Phi_d, X_d, D_d, para_d); cudaThreadSynchronize ();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        // Anzahl der Verletzten neu bestimmen und auf Host speichern
        sumUp<<<1,1>>> (Injured_d,N, NInjured_d); cudaThreadSynchronize();  error = cudaGetLastError ();
        error = cudaMemcpy(&NInjured, NInjured_d, sizeof(int), cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        // Summe der Kräfte berechnen
        sqrt_fact = sqrt(tstep/DefaultDeltaT);
        sumForces<<<dimGrid, dimBlock >>> (fsumx_d, fsumy_d, tStepVector_d, sqrt_fact, VX_d, VY_d, V0of_d, Phi_d, fpairx_d, fwallx_d, fwpointx_d, fpairy_d, fwally_d, fwpointy_d, fsmokex_d, fsmokey_d, fcolx_d, fcoly_d,N, para_d);
        cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        // neuen timeStep bestimmen und auf Host speichern
        getMinTimeStep <<<1,1>>> (tStepVector_d,N, timeStep_d);  cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        error = cudaMemcpy(&tstep, timeStep_d, sizeof(float), cudaMemcpyDeviceToHost); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        // LOG(WARNING) << "newTimeStep :" << tstep;


        NewVelocity <<<dimGrid, dimBlock >>> (vxnew_d,  vynew_d,   fsumx_d,   fsumy_d,   VX_d,   VY_d,   Injured_d, N, tStepVector_d , para_d); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        getNewValues <<<dimGrid, dimBlock >>> (Xprev_d, X_d, Yprev_d, Y_d, VY_d, VX_d, NinRoomVektor_d, tstep, N, para_d); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);

        // benötigte Werte vorher runterkopieren
        // calcNewBookKeepingOnHost ()
        // veränderte Werte wieder hochkopieren

        // storeNewVelocity
        storeNewVelocity <<<dimGrid, dimBlock >>> (VX_d, VY_d, V_d, Vdir_d,  Phi_d, X_d, Y_d, D_d, W_d,  vxnew_d, vynew_d, para_d, N, YS); cudaThreadSynchronize();
        error = cudaGetLastError(); CHECK_EQ(cudaSuccess, error) << "Error: " << cudaGetErrorString(error);
        LOG (INFO) << "Updateschritt beendet." ;
        
		// Ende des Update Schrittes
        SimTime[UpdNum+1] = SimTime[UpdNum] + tstep;
        UpdNum++;
        
    } while(( UpdNum < MaxUpdNum ) &&( SimTime[UpdNum] < MaxSimTime ));
    LOG (WARNING) << "Schleife beendet.";
    
    return 0;
}