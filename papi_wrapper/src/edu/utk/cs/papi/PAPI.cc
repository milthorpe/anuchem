
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <papi.h>
#include "PAPI.h"

namespace edu {
    namespace utk {
        namespace cs {
            namespace papi {
                PAPI* PAPI::_make() {
                    return new PAPI();
                }

                PAPI::PAPI() {
                    int num_hwcntrs;
                    if ((num_hwcntrs = PAPI_num_counters()) <= PAPI_OK)  
                        printerror(__FILE__, __LINE__, "PAPI_num_counters", num_hwcntrs);
                    values = new long long[num_hwcntrs];
                }

                void PAPI::initialize() {
                    int retval = PAPI_library_init(PAPI_VER_CURRENT);

                    if (retval != PAPI_VER_CURRENT && retval > 0) {
                        fprintf(stderr,"PAPI library version mismatch!\n");
                        exit(1);
                    }

                    if (retval < 0) {
                        fprintf(stderr, "Initialization error!\n");
                        exit(1);
                    }
                }

                void PAPI::shutDown() {
                    if (EventSet != NULL) destroyEventSet();
                    PAPI_shutdown();
                }

                void PAPI::startFlops() {
                    createEventSet();
                    addEvent(PAPI_TOT_INS);
                    addEvent(PAPI_FP_OPS);
                    start();
                }

                void PAPI::start() {
                    int retval;
                    /* Start counting events in the Event Set */
                    if (retval=PAPI_start(EventSet) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_start", retval);
                }

                void PAPI::stop() {
                    int retval;
                    /* Read the counting events in the Event Set */
                    if (retval=PAPI_stop(EventSet, values) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_stop", retval);
                }

                void PAPI::createEventSet() {
                    int retval;
                    if (retval=PAPI_create_eventset(&EventSet) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_create_eventset", retval);
                }

                void PAPI::addEvent(int EventCode) {
                    int retval;
                    if (retval=PAPI_add_event(EventSet, EventCode) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_add_event", retval);
                }

                void PAPI::destroyEventSet() {
                    int retval;
                    if (retval=PAPI_cleanup_eventset(EventSet) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_cleanup_eventset", retval);
                    if (retval=PAPI_destroy_eventset(&EventSet) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_destroy_eventset", retval);
                }

                int64_t PAPI::getCounter(int i) {
                    return values[i];
                }
            }
        }
    }
}

static void printerror(char *file, int line, char *call, int retval) {
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

