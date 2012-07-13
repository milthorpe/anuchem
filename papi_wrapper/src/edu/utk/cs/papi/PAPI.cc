
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
                    totals = new long long[num_hwcntrs];
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
                    addEvent(PAPI_FP_OPS);
                    resetTotals();
                }

                void PAPI::printFlops() {
                    printf("cycles: %16lld total ins: %16lld FP ins: %16lld FLOPS: %16lld\n", cycles, totals[0], totals[1], totals[2]);
                }

                void PAPI::countMemoryOps() {
                    if (EventSet != PAPI_NULL) destroyEventSet();
                    createEventSet();
                    addEvent(PAPI_TOT_INS);
                    addEvent(PAPI_LD_INS);
                    addEvent(PAPI_SR_INS);
                    //addEvent(PAPI_L1_LDM);
                    //addEvent(PAPI_L1_STM);
                    resetTotals();
                }

                void PAPI::printMemoryOps() {
                    printf("cycles: %16lld total ins: %16lld LD ins: %16lld SR ins: %16lld\n", cycles, totals[0], totals[1], totals[2]);
                    //printf("cycles: %16lld total ins: %16lld LD ins: %16lld SR ins: %16lld L1 LD miss: %16lld L1 ST miss: %16lld\n", cycles, totals[0], totals[1], totals[2], totals[3], totals[4]);
                }

                void PAPI::start() {
                    cycles -= PAPI_get_real_cyc();
                    int retval;
                    if ((retval=PAPI_start(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_start", retval);
                }

                void PAPI::stop() {
                    cycles += PAPI_get_real_cyc();
                    int retval;
                    if ((retval=PAPI_stop(EventSet, values)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_stop", retval);
                    for (int i=0; i<num_hwcntrs; i++) {
                        totals[i] += values[i];
                    }
                }

                void PAPI::reset() {
                    resetTotals();
                    int retval;
                    /* Reset the counting events in the Event Set */
                    if ((retval=PAPI_reset(EventSet)) != PAPI_OK)
                        printerror(__FILE__, __LINE__, "PAPI_reset", retval);

                }

                int64_t PAPI::getCounter(int i) {
                    return totals[i];
                }

                void PAPI::resetCounter(int i) {
                    totals[i] = 0;
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

                void PAPI::resetTotals() {
                    cycles = 0LL;
                    bzero(totals, num_hwcntrs*sizeof(double));
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

