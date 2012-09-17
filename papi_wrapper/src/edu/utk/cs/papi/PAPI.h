

namespace edu {
    namespace utk {
        namespace cs {
            namespace papi {
                class PAPI {
                    public:
                        static PAPI* _make();
                        void initialize();
                        void countFlops();
                        void printFlops();
                        void countMemoryOps();
                        void printMemoryOps();
                        void start();
                        void stop();
                        void reset();
                        int64_t getCounter(int i);
                        void resetCounter(int i);
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


