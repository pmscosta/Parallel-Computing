#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <chrono>
#include <cstdint>
#include <omp.h>
#include <fstream>
#include <vector>
#include <map>

using namespace std;

ofstream timeLogger;
int START_EXP = 25;
int END_EXP = 32;

int sieve_seq(u_long lastNumber)
{
    const u_long lastNumberSqrt = (int)sqrt((double)lastNumber);

    u_long memorySize = (lastNumber - 1) / 2;

    // 0 -> it's prime, 1 it's not prime
    bool *isPrime = new bool[memorySize + 1];

    for (u_long i = 3; i <= lastNumberSqrt; i += 2)
    {
        if (!isPrime[i / 2])
            for (u_long j = i * i; j <= lastNumber; j += 2 * i)
                isPrime[j >> 1] = 1;
    }

    int found = lastNumber >= 2 ? 1 : 0;

    for (u_long i = 1; i <= memorySize; i++)
        found += !isPrime[i];

    delete[] isPrime;
    return found;
}

int sieve_omp(u_long lastNumber)
{

    omp_set_num_threads(omp_get_num_procs());

    const u_long lastNumberSqrt = (int)sqrt((double)lastNumber);

    u_long memorySize = (lastNumber - 1) / 2;

    // 0 -> it's prime, 1 it's not prime
    bool *isPrime = new bool[memorySize + 1];

#pragma omp parallel for schedule(dynamic)
    for (u_long i = 3; i <= lastNumberSqrt; i += 2)
    {
        if (!isPrime[i / 2])
            for (u_long j = i * i; j <= lastNumber; j += 2 * i)
                isPrime[j >> 1] = 1;
    }

    int found = lastNumber >= 2 ? 1 : 0;

#pragma omp parallel for reduction(+ \
                                   : found)
    for (u_long i = 1; i <= memorySize; i++)
        found += !isPrime[i];

    delete[] isPrime;
    return found;
}

int average_counter(int numTimes, int f(u_long))
{
    u_long start = pow(2, START_EXP);
    u_long end = pow(2, END_EXP);

    // 3 iterations. don't know how many are needed
    u_long step = (end - start) / 4;

    map<u_long, double> times;

    for (int i = 0; i < numTimes; i++)
    {
        for (u_long val = start; val <= end; val += step)
        {
            auto t1 = std::chrono::high_resolution_clock::now();
            int found = f(val);
            auto t2 = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

            double seconds = duration / 1000000.0;

            // cout << "Found: " << found << " primes in " << seconds << " seconds." << endl;

            times[val] += seconds;
        }
    }

    double average = 0;

    for (auto const &entry : times)
    {
        cout << "For " << entry.first << ", average: " << entry.second / 4.0 << " seconds\n";
        timeLogger << entry.first << "," << entry.second / 4.0 << endl;
    }
}

int main(int argc, char *argv[])
{
    timeLogger.open("timeLogger.txt", std::ios_base::app);
    int counter = atoi(argv[1]);
    average_counter(counter, sieve_seq);
    timeLogger.close();
    return 0;
}
