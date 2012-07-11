

namespace edu {
    namespace utk {
        namespace cs {
            namespace papi {
                class PAPI {
                    public:
                        static PAPI* _make();
                        PAPI();
                        ~PAPI();
                        void initialize();
                        void countFlops();
                        void printFlops();
                        void countMemoryOps();
                        void printMemoryOps();
                        void startCount();
                        void stopCount();
                        void resetCount();
                        int64_t getCounter(int i);
                        void createEventSet();
                        void addEvent(int EventCode);
                        void destroyEventSet();
                        void shutDown();
				    private:
                        int EventSet;
                        int num_hwcntrs;
                        long long cycles;
					    long long *values;
					    long long *totals;
                        void resetTotals();
                };
            }
        }
    }
}


