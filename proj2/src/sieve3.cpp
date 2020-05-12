#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(((id) + 1), p, n) - 1)


void usage(void)
{
    std::cout << "sieve <max number>" << std::endl;
    std::cout << "<max number> range between 2 and N." << std::endl;
}

ulong singleBlock(ulong from, ulong to);

int main(int argc, char *argv[])
{
    double elapsed_time;

    MPI_Init(&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    int process_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_id);

    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    if (argc != 2)
    {
        if (process_id == 0)
        {
            usage();
            MPI_Finalize();
            exit(1);
        }
    }

    ulong range_max = atol(argv[1]);

    // ulong sqrtn = ceil(sqrt((ulong)range_max));

    ulong low_value = 2 + BLOCK_LOW(process_id, num_processes, range_max -1);
    ulong high_value = 2+ BLOCK_HIGH(process_id, num_processes, range_max-1);

    // printf("%ld %ld\n", low_value, high_value);

    ulong from;
    ulong found = 0;
    ulong count = 0;
    const ulong slice = 256000;

    for (from = low_value; from <= high_value; from += slice)
    {
        ulong to = from + slice;

        if (to > high_value)
            to = high_value;

        found += singleBlock(from, to);
    
    }

    MPI_Reduce(&found, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    elapsed_time += MPI_Wtime();

    if (process_id == 0)
    {
        printf("Found %lu primes\n", count);
        printf("Time: %3.3f seconds \n", elapsed_time);

    }

    MPI_Finalize();

    return 0;
}

ulong singleBlock(ulong from, ulong to)
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

        ulong  minJ = ((from + i - 1) / i) * i;
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