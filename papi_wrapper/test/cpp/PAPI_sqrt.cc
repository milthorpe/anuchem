/*
 * A simple test of sqrt / div FP instructions.
 */
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <memory.h>
#include <malloc.h>
#include <math.h>
#include "papi.h"

#define ITERS 500000


static void printerror(const char *file, int line, const char *call, int retval);

int main(int argc, char **argv) {
    int EventSet=PAPI_NULL;
    long long *values = new long long[2];
    int retval = PAPI_library_init(PAPI_VER_CURRENT);

    if ((retval=PAPI_create_eventset(&EventSet)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_create_eventset", retval);
    if ((retval=PAPI_add_event(EventSet, PAPI_TOT_INS)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_add_event", retval);
    if ((retval=PAPI_add_event(EventSet, PAPI_FP_INS)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_add_event", retval);

    long cycles = -PAPI_get_real_cyc();
    if ((retval=PAPI_start(EventSet)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_start", retval);

    double j = 2.0;
    for (int i=0; i<ITERS; i++) {
        j += sqrt(j); // j += 1.0 / j;
    }

    cycles += PAPI_get_real_cyc();
    if ((retval=PAPI_stop(EventSet, values)) != PAPI_OK)
        printerror(__FILE__, __LINE__, "PAPI_stop", retval);
                    printf("cycles: %16ld total ins: %16lld FP ins: %16lld\n", cycles, values[0], values[1]);

    printf("total = %f\n", j);

    PAPI_shutdown();
    exit(0);
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

