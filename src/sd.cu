#include "sd_lib.cu"
#include "update.cu"



int main( int NArg, char * ArgStr[] )
{
	google::InitGoogleLogging(ArgStr[0]);
    Start_Bare( NArg, ArgStr ); // Ursprünglich war hier Start_Demo (aber da wurde nur Start_Bare aufgerufen
    Init_Demo();
	
	parameter para_h;
	prepareParameter (&para_h);

	LOG (WARNING) << "Hauptschleife wird gestartet.";
    do {
        Pic();
        Save_Demo();
        Upd(para_h);
        // ReInit habe ich entfernt, weil im Projekt da kein Wert darauf gelegt wird

    } while( /* (N>0)&& */ ( UpdNum < MaxUpdNum ) &&( SimTime[UpdNum] < MaxSimTime ));
	LOG (WARNING) << "Hauptschleife beendet mit " << MaxUpdNum << " Update Schritten";
    Shutdown_Demo();

    return 0;
}
