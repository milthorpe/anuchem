
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <papi.h>
#include "PAPI.h"

static void printerror(const char *file, int line, const char *call, int retval);

namespace edu {
    namespace utk {
        namespace cs {
            namespace papi {
                PAPI* PAPI::_make() {
                    return new PAPI();
                }

                PAPI::PAPI() {
                    EventSet=PAPI_NULL;
                    if ((num_hwcntrs = PAPI_num_counters()) <= PAPI_OK)  
                        printerror(__FILE__, __LINE__, "PAPI_num_counters", num_hwcntrs);
                    values = new long long[num_hwcntrs];
                    total = new long long[num_hwcntrs];
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

                PAPI::~PAPI() {
                    shutDown();
                }

                void PAPI::countFlops() {
                    if (EventSet != PAPI_NULL) destroyEventSet();
                    createEventSet();
                    addEvent(PAPI_TOT_INS);
                    addEvent(PAPI_FP_INS);
                    //addEvent(PAPI_FP_OPS);
                    startCount();
                }

                void PAPI::printFlops() {
                    printf("total cyc: %16lld total ins: %16lld total FP ins: %16lld total FLOPS: %16lld\n", cycles, values[0], values[1], values[2]);
                }

                void PAPI::countMemoryOps() {
                    if (EventSet != PAPI_NULL) destroyEventSet();
                    createEventSet();
                    addEvent(PAPI_TOT_INS);
                    addEvent(PAPI_LD_INS);
                    //addEvent(PAPI_SR_INS);
                    startCount();
                }

                void PAPI::printMemoryOps() {
                    printf("total cyc: %16lld total ins: %16lld total LD ins: %16lld total SR ins: %16lld\n", cycles, values[0], values[1], values[2]);
                }

                void PAPI::startCount() {
                    int retval;
                    if ((retval=PAPI_start(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_start", retval);
                }

                void PAPI::stopCount() {
                    cycles += PAPI_get_real_cyc();
                    int retval;
                    if ((retval=PAPI_accum(EventSet, values)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_accum", retval);
                }

                void PAPI::resumeCount() {
                    cycles -= PAPI_get_real_cyc();
                    int retval;
                    /* Reset the counting events in the Event Set */
                    if ((retval=PAPI_reset(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_reset", retval);
                }

                void PAPI::resetCount() {
                    cycles = -PAPI_get_real_cyc();
                    bzero(values, num_hwcntrs*sizeof(double));
                    int retval;
                    /* Reset the counting events in the Event Set */
                    if ((retval=PAPI_reset(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_reset", retval);

                }

                int64_t PAPI::getCounter(int i) {
                    return values[i];
                }

                void PAPI::createEventSet() {
                    int retval;
                    if ((retval=PAPI_create_eventset(&EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_create_eventset", retval);
                }

                void PAPI::addEvent(int EventCode) {
                    int retval;
                    if ((retval=PAPI_add_event(EventSet, EventCode)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_add_event", retval);
                }

                void PAPI::destroyEventSet() {
                    int retval;
                    if ((retval=PAPI_cleanup_eventset(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_cleanup_eventset", retval);
                    if ((retval=PAPI_destroy_eventset(&EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_destroy_eventset", retval);
                }

                void PAPI::shutDown() {
                    if (EventSet != PAPI_NULL) destroyEventSet();
                    PAPI_shutdown();
                }
            }
        }
    }
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

