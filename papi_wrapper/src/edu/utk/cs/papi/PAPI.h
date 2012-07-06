

namespace edu {
    namespace utk {
        namespace cs {
            namespace papi {
                class PAPI {
                    public:
                        PAPI();
                        static PAPI* _make();
                        void initialize();
                        void shutDown();
                        void startFlops();
                        void start();
                        void stop();
                        void createEventSet();
                        void addEvent(int EventCode);
                        void destroyEventSet();
                        int64_t getCounter(int i);
				    private:
                        int EventSet;
					    long long *values;
                };
            }
        }
    }
}

static void printerror(char *file, int line, char *call, int retval);

