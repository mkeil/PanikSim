#ifndef ABS
#define	ABS(a)	        ( (a) >= 0.0 ? (a) : (-1.0)*(a) )
#endif

#ifndef SIGN
#define SIGN(a)        ( (a) >= 0.0 ? (1.0) : (-1.0) )
#endif

#ifndef SQR
#define SQR(a)          ((a)*(a))
#endif

#ifndef MAX
#define	MAX(a, b)	( (a) >= (b) ? (a) : (b) )
#endif

#ifndef MIN
#define	MIN(a, b)	( (a) <= (b) ? (a) : (b) )
#endif

#ifndef LIM
#define	LIM(L1,a,L2)	( (L1) <= (L2) ? ( MIN( (MAX((L1),(a))), (L2) )) : ( MIN( (MAX((L2),(a))), (L1) )) )
#endif

#ifndef PI
#define PI 3.14159265358979
#endif

#ifndef _E
#define _E(s) {fprintf(stderr,"%s",s);fflush(stderr);}
#endif

/*#ifndef _ED
 *char _MY_STRING_FOR_ED_[1000];
 *#define _ED(s) {fprintf(stderr,"%s\t");fflush(stderr);sprintf(_MY_STRING_FOR_ED_,"date");system(_MY_STRING_FOR_ED_);}
 *#endif
 */

#ifndef _O
#define _O(s) {fprintf(stdout,"%s",s);fflush(stdout);}
#endif
