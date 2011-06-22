// foo
/* self-driven / panic / main file for demo version
   MainSwitch options:
   0 simple demo --> output on X11
   1 eps images instead of X11
   2 creating data file for java
   */

/*********************/


#include "sd_lib.cu"
#include "update.cu"

int main( int NArg, char * ArgStr[] )
{

    Start_Bare( NArg, ArgStr ); // Ursprünglich war hier Start_Demo (aber da wurde nur Start_Bare aufgerufen
    Init_Demo();




    do {
        Pic();
        Save_Demo();
        Upd();
        ReInit();

    } while( /* (N>0)&& */ ( UpdNum < MaxUpdNum ) &&( SimTime[UpdNum] < MaxSimTime ));

    Shutdown_Demo();




    return 0;
}
