#define _GLIBCXX_PARALLEL
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream>
#include <stdio.h>
#include <omp.h>

using namespace std;

static long microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}


void setToOne(int& i, int val)
{
  i = val;
}

void coutVector(const std::vector<int>& myVector)
{
  const int size = myVector.size();
  for (int i=0; i!=size; ++i)
  {
    std::cout << i << " : " << myVector[i] << std::endl;
  }
}

const int ITERS = 100000;

int main() {
    std::vector<int> myVector(8000);

    long start = microTime();

    for (int i=0; i<ITERS; i++) {
        for (int j=0; j<8000; j++) {
            myVector[j] = j;
        }
    }

    long stop = microTime();

    printf("serial for loop %.3f ms\n", (stop-start) / 1.0e3 / ITERS);

    start = microTime();

    for (int i=0; i<ITERS; i++) {
        std::for_each(myVector.begin(), myVector.end(), [](int &n){ n++; });
    }

    stop = microTime();

    printf("std::for_each %.3f ms\n", (stop-start) / 1.0e3 / ITERS);

    //coutVector(myVector);
    return 0;
}

