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

u_long sieve_omp_slices(u_long lastNumber, u_long slices)
{

    omp_set_num_threads(omp_get_num_procs());

    const u_long lastNumberSqrt = (u_long)sqrt((u_long)lastNumber);

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

    u_long found = lastNumber >= 2 ? 1 : 0;

#pragma omp parallel for reduction(+ \
                                   : found)
    for (u_long i = 1; i <= memorySize; i++)
        found += !isPrime[i];

    delete[] isPrime;
    return found;
}

u_long sieve_omp_odd_single_block(const u_long from, const u_long to)
{
    const ulong memorySize = (to - from + 1) / 2.0;

    char *isPrime = new char[memorySize];

    for (ulong i = 0; i < memorySize; i++)
        isPrime[i] = 1;

    for (ulong i = 3; i * i <= to; i += 2)
    {
        if (i >= 3 * 3 && i % 3 == 0)
            continue;
        if (i >= 5 * 5 && i % 5 == 0)
            continue;
        if (i >= 7 * 7 && i % 7 == 0)
            continue;
        if (i >= 11 * 11 && i % 11 == 0)
            continue;
        if (i >= 13 * 13 && i % 13 == 0)
            continue;

        ulong minJ = ((from + i - 1) / i) * i;
        if (minJ < i * i)
            minJ = i * i;
        // start value must be odd
        if ((minJ & 1) == 0)
            minJ += i;

        // find all odd non-primes
        for (ulong j = minJ; j <= to; j += 2 * i)
        {
            ulong index = j - from;
            isPrime[index / 2] = 0;
        }
    }

    ulong found = 0;
    for (ulong i = 0; i < memorySize; i++)
        found += isPrime[i];
    // 2 is not odd => include on demand
    if (from <= 2)
        found++;

    return found;
}

u_long sieve_omp_blockwise_parallel(u_long lastNumber, u_long sliceSize)
{
    omp_set_num_threads(omp_get_num_procs());

    u_long found = 0;

#pragma omp parallel for reduction(+ \
                                   : found)
    for (u_long from = 2; from <= lastNumber; from += sliceSize)
    {
        u_long to = from + sliceSize;
        if (to > lastNumber)
        {
            to = lastNumber;
        }

        found += sieve_omp_odd_single_block(from, to);
    }

    return found;
}

u_long sieve_omp_blockwise(u_long lastNumber, u_long sliceSize)
{
    u_long found = 0;

    for (u_long from = 2; from <= lastNumber; from += sliceSize)
    {
        u_long to = from + sliceSize;
        if (to > lastNumber)
        {
            to = lastNumber;
        }

        found += sieve_omp_odd_single_block(from, to);
    }

    return found;
}

int average_counter(int numTimes, int f(u_long))
{
    u_long start = pow(2, START_EXP);
    u_long end = pow(2, END_EXP);

    // 3 iterations. don't know how many are needed
    u_long step = (end - start) / 4;

    map<u_long, vector<double>> times;

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

            times[val].push_back(seconds);
        }
    }

    double average = 0;

    for (auto const &entry : times)
    {
        double total = 0;
        for (double time : entry.second)
        {
            total += time;
        }

        double avg_time = (double)(total / entry.second.size());
        cout << "For " << entry.first << ", average: " << avg_time << " seconds\n";
        timeLogger << entry.first << "," << avg_time << endl;
    }
}

int average_counter_slices(int numTimes, u_long f(u_long, u_long))
{
    u_long start = pow(2, START_EXP);
    u_long end = pow(2, END_EXP);

    // 3 iterations. don't know how many are needed
    u_long step = (end - start) / 4;

    map<u_long, vector<double>> times;
    ulong slice = 1024 * 512;
    for (int i = 0; i < numTimes; i++)
    {
        for (u_long val = start; val <= end; val += step)
        {

            auto t1 = std::chrono::high_resolution_clock::now();
            ulong found = f(val, slice);
            auto t2 = std::chrono::high_resolution_clock::now();

            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

            double seconds = duration / 1000000.0;

            cout << "Found: " << found << " primes in " << seconds << " seconds." << endl;

            times[val].push_back(seconds);
        }
    }

    double average = 0;

    for (auto const &entry : times)
    {
        double total = 0;
        for (double time : entry.second)
        {
            total += time;
        }

        double avg_time = (double)(total / entry.second.size());

        cout << "For " << entry.first << " slice " << slice << ", average: " << avg_time << " seconds\n";
        timeLogger << entry.first << "," << avg_time << endl;
    }
}

int main(int argc, char *argv[])
{
    int counter;
    // timeLogger.open("timeLogger-seq.txt", std::ios_base::app);
    // counter = atoi(argv[1]);
    // average_counter(counter, sieve_seq);
    // timeLogger.close();

    // timeLogger.open("timeLogger-omp.txt", std::ios_base::app);
    // counter = atoi(argv[1]);
    // average_counter(counter, sieve_omp);
    // timeLogger.close();

    // timeLogger.open("timeLogger-omp-blockwise.txt", std::ios_base::app);
    // counter = atoi(argv[1]);
    // average_counter_slices(counter, sieve_omp_blockwise);
    // timeLogger.close();

    timeLogger.open("timeLogger-omp-blockwise-parallel.txt", std::ios_base::app);
    counter = atoi(argv[1]);
    average_counter_slices(counter, sieve_omp_blockwise_parallel);
    timeLogger.close();

    return 0;
}
