#define _GLIBCXX_PARALLEL
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include <sys/time.h>
#include <iostream>
#include <stdio.h>

using namespace std;

static long microTime() {
    struct ::timeval tv;
    gettimeofday(&tv, NULL);
    return (long)(tv.tv_sec * 1000000LL + tv.tv_usec);
}

void setToOne(int& i)
{
  i = 1;
}

void coutVector(const std::vector<int>& myVector)
{
  const int size = myVector.size();
  for (int i=0; i!=size; ++i)
  {
    std::cout << i << " : " << myVector[i] << std::endl;
  }
}


int main() {
    std::vector<int> myVector(8000);

    long start = microTime();

    for (int i=0; i<10000; i++) {
        for (int j=0; j<8000; j++) {
            myVector[j] = 1;
        }
    }

    long stop = microTime();

    printf("serial for loop %.3f ms\n", (stop-start) / 1.0e4);

    start = microTime();

    for (int i=0; i<10000; i++) {
        std::for_each(myVector.begin(),myVector.end(), setToOne);
    }

    stop = microTime();

    printf("std::for_each %.3f ms\n", (stop-start) / 1.0e4);

    //coutVector(myVector);
    return 0;
}

