/* readpar_v_2000_03_02.c

   reading parameter file with this format:
   1st column     #
   2nd            one of these characters: i f s (integer/float/string)
   3rd            one of these characters: 0 1 (1:interactive parameter,0:not)
   4th            parameter name
   5th            parameter value
*/

#include <stdio.h>
#include <string.h>
#define READPAR_EXIT {fprintf(stderr,"readpar EXIT\n");fflush(stderr);exit(-1);}


void readpar ( char *sw, char *ifn,
               int *intValue[], char *intName[], int intNum,
               float *floatValue[], char *floatName[], int floatNum,
               char *stringValue[], char *stringName[], int stringNum )
{

    /* sw: switch = "start" or "re"
       ifn: input file name
       each array starts with the 0. element
       */


    FILE *ifp;
    int ii,i,*intFound,*floatFound,*stringFound,exitFlag,
        interactive,tmpInt;
    char tmpString[100], tmpStringVal[100],c;
    float tmpFloat;


    /* malloc, init, etc */
    intFound = ivector(0,intNum-1);
    floatFound = ivector(0,floatNum-1);
    stringFound = ivector(0,stringNum-1);
    for(i=0; i<intNum; i++) {
        intFound[i]=0;
    }
    for(i=0; i<floatNum; i++) {
        floatFound[i]=0;
    }
    for(i=0; i<stringNum; i++) {
        stringFound[i]=0;
    }

    if( !(ifp = fopen(ifn, "r")) ) {
        fprintf(stderr,"readpar: Couldn't open \"%s\" for reading\n",ifn);
        READPAR_EXIT;
    }


    for(ii=0; ii<intNum+floatNum+stringNum; ii++) {
        while(getc(ifp)!=0x23) {}; /* 0x23 = '#' */

        /* type of parameter */
        fscanf(ifp,"%s",&c);
        switch(c) {
        case 0x69: { /* 0x69 = 'i' */
                exitFlag=0;
                fscanf(ifp,"%d %s %d",&interactive,tmpString,&tmpInt);
                for(i=0; (exitFlag==0)&&(i<intNum); i++) {
                    if(strcmp(intName[i],tmpString)==0) {
                        intFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            *(intValue[i])=tmpInt;
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this integer parameter: %s\n",tmpString);
                }
                break;
            }
        case 0x66: { /* 0x66 = 'f' */
                exitFlag=0;
                fscanf(ifp,"%d %s %f",&interactive,tmpString,&tmpFloat);
                for(i=0; (exitFlag==0)&&(i<floatNum); i++) {
                    if(strcmp(floatName[i],tmpString)==0) {
                        floatFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            *(floatValue[i])=tmpFloat;
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this float parameter: %s\n",tmpString);
                }
                break;
            }
        case 0x73: { /* 0x73 = 's' */
                exitFlag=0;
                fscanf(ifp,"%d %s %s",&interactive,tmpString,tmpStringVal);
                for(i=0; (exitFlag==0)&&(i<stringNum); i++) {
                    if(strcmp(stringName[i],tmpString)==0) {
                        stringFound[i]=1;
                        exitFlag=1;
                        if(  (strcmp(sw,"start")==0)
                                ||((strcmp(sw,"re")==0)&&(interactive==1))
                          ) {
                            strcpy(stringValue[i],tmpStringVal);
                        }
                    }
                }
                if(exitFlag==0) {
                    fprintf(stderr,"readpar WARNING: don't need this string parameter: %s\n",tmpString);
                }
                break;
            }
        }
    }





    /* checking whether all parameters have been found */
    for(i=0; i<intNum; i++) {
        if(intFound[i]==0) {
            fprintf(stderr,"readpar ERROR: integer parameter %s not found in %s\n",intName[i],ifn);
            READPAR_EXIT;
        }
    }
    for(i=0; i<floatNum; i++) {
        if(floatFound[i]==0) {
            fprintf(stderr,"readpar ERROR: float parameter %s not found in %s\n",floatName[i],ifn);
            READPAR_EXIT;
        }
    }
    for(i=0; i<stringNum; i++) {
        if(stringFound[i]==0) {
            fprintf(stderr,"readpar ERROR: stringparameter %s not found in %s\n",stringName[i],ifn);
            READPAR_EXIT;
        }
    }

    free_ivector(intFound,0,intNum-1);
    free_ivector(floatFound,0,floatNum-1);
    free_ivector(stringFound,0,stringNum-1);
}
