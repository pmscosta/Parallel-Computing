#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <chrono>
#include <cstdint>
#include <omp.h>

using namespace std;

int sieve(uint32_t lastNumber, bool useOpenMP);

int main()
{

    auto t1 = std::chrono::high_resolution_clock::now();
    int found = sieve(1 * 1000 * 1000 * 1000LL, true);
    auto t2 = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    double seconds = duration / 1000000.0;

    cout << "Found: " << found << " primes in " << seconds << " seconds." << endl;

    return 0;
}

int sieve(uint32_t lastNumber, bool useOpenMP)
{

    omp_set_num_threads(useOpenMP ? omp_get_num_procs() : 1);

    const uint32_t lastNumberSqrt = (int)sqrt((double)lastNumber);

    uint32_t memorySize = (lastNumber - 1) / 2;

    // 0 -> it's prime, 1 it's not prime
    char *isPrime = new char[memorySize + 1];

    #pragma omp parallel for schedule(dynamic)
        for (uint32_t i = 3; i <= lastNumberSqrt; i += 2){
            if (!isPrime[i / 2])
                for (uint32_t j = i * i; j <= lastNumber; j += 2 * i)
                    isPrime[j / 2] = 1;
        }

    int found = lastNumber >= 2 ? 1 : 0;

    #pragma omp parallel for reduction(+ : found)
        for (uint32_t i = 1; i <= memorySize; i++)
            found += !isPrime[i];

    delete[] isPrime;
    return found;
}