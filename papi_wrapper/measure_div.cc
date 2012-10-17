#include <stdlib.h>
#include <sys/time.h> 
#include <math.h>
#include <stdio.h>
#include <papi.h>
#include <string.h>
#include <strings.h>

#define USEPAPI 1

using namespace std;

static long microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}

static void printerror(const char *file, int line, const char *call, int retval);

static double doAMillionDivs(double a[]) {
#ifdef USEPAPI
    int EventSet=PAPI_NULL;
    int retval = PAPI_library_init(PAPI_VER_CURRENT);

    if (retval != PAPI_VER_CURRENT && retval > 0) {
        fprintf(stderr,"PAPI library version mismatch!\n");
        exit(1);
    }

    if (retval < 0) {
        fprintf(stderr, "Initialization error!\n");
        exit(1);
    }
    PAPI_create_eventset(&EventSet);
    PAPI_add_event(EventSet, PAPI_TOT_INS);
    PAPI_add_event(EventSet, PAPI_FP_INS);
    PAPI_add_event(EventSet, PAPI_FDV_INS);
    long long *totals = new long long[3];
    bzero(totals, 3*sizeof(long long));
    long long cycles = -PAPI_get_real_cyc();
    if ((retval=PAPI_start(EventSet)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_start", retval);
#endif

    for (int i=0; i<1000000; i++) {
        a[i] = a[i] / 1.5;
    }

#ifdef USEPAPI
    cycles += PAPI_get_real_cyc();
    PAPI_stop(EventSet, totals);
    printf("cycles: %16lld total ins: %16lld FP ins: %16lld div: %16lld\n", cycles, totals[0], totals[1], totals[2]);
    printf("cycles/div %f\n",  cycles / 1e6);
#endif

}

static void printerror(const char *file, int line, const char *call, int retval) {
    printf("%s\tFAILED\nLine # %d\n", file, line);
    if ( retval == PAPI_ESYS ) {
        char buf[128];
        memset( buf, '\0', sizeof(buf) );
        sprintf(buf, "System error in %s:", call );
        perror(buf);
    }
    else if ( retval > 0 ) {
        printf("Error calculating: %s retval %d\n", call, retval );
    }
    else {
        char errstring[PAPI_MAX_STR_LEN];
        PAPI_perror(retval, errstring, PAPI_MAX_STR_LEN );
        printf("Error in %s: %s\n", call, errstring );
    }
    printf("\n");
    exit(1);
}

int main(int argc, char* argv[]) {
    double a[1000000];
    for (int i=0; i<1000000; i++) {
        a[i] = i+1.0;
    }
    long start = microTime();
    doAMillionDivs(a);
    long stop = microTime();
    printf("%.3g s\n", (stop-start) / 1.0e6);

    double res=0.0;
    for (int i=0; i<1000000; i++) {
        res += a[i];
    }
    printf("%f\n", res);

    return 0;
}
