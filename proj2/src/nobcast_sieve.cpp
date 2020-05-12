#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <vector>

#define BLOCKLOW(id, p, n) ((id) * (n) / (p))
#define BLOCKHIGH(id, p, n) (BLOCKLOW(((id) + 1), p, n) - 1)
#define BLOCKSIZE(id, p, n) ((BLOCKLOW(((id) + 1), p, n)) - (BLOCKLOW(id, p, n)))

int main(int argc, char *argv[])
{

    int id;
    int p; // number of processes
    unsigned long n;
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    double time = MPI_Wtime();

    n = atol(argv[1]);

    unsigned long sqrt_n = ceil(sqrt((double)n));

    std::vector<bool> pre_marked(sqrt_n + 1, false);

    pre_marked.at(0) = true;
    pre_marked.at(1) = true;

    int pre_k = 2;

    // marcar todos os even dado que o 2 é o único even que é primo
    do
    {
        int base = pre_k * pre_k;
        for (int i = base; i <= sqrt_n; i += pre_k)
            pre_marked[i] = true;

        while (pre_marked[++pre_k])
        {
        };

    } while (pre_k * pre_k <= sqrt_n);

    std::vector<int> kset;
    for (int i = 3; i <= sqrt_n; ++i)
    {
        if (!pre_marked[i])
        {
            kset.push_back(i);
        }
    }

    unsigned long low_value = 2 + BLOCKLOW(id, p, n - 1);
    unsigned long high_value = 2 + BLOCKHIGH(id, p, n - 1);
    unsigned long block_size = BLOCKSIZE(id, p, n - 1);

    if (low_value % 2 == 0)
    {
        if (high_value % 2 == 0)
        {
            block_size = (unsigned long)floor((unsigned long)block_size / 2.0);
            high_value--;
        }
        else
        {
            block_size = block_size / 2;
        }
        low_value++;
    }
    else
    {
        if (high_value % 2 == 0)
        {
            block_size = block_size / 2;
            high_value--;
        }
        else
        {
            block_size = (unsigned long)ceil((unsigned long)block_size / 2.0);
        }
    }

    std::vector<bool> marked(block_size, false);

    unsigned long first_index;

    if (id == 0)
    {
        first_index = 0;
    }

    unsigned long kindex = 0;

    int k = kset[kindex];

    unsigned long count = 1;

    do{
        if (k >= low_value)
        {
            first_index = ((k - low_value) / 2) + k;
        }
        else if (k * k > low_value)
        {
            first_index = (k * k - low_value) / 2;
        }
        else
        {
            if (low_value % k == 0)
            {
                first_index = 0;
            }
            else
            {
                first_index = 1;
                while ((low_value + (2 * first_index)) % k != 0)
                    ++first_index;
            }
        }

        for (int i = first_index; i < block_size; i += k)
        {
            marked[i] = true;
        }

        k = kset[++kindex];
    } while (k * k <= n && kindex < kset.size());

    int local_count = 0;

    for (int i = 0; i < block_size; i++)
    {
        if (!marked[i])
            local_count++;
    }

    MPI_Reduce(&local_count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (id == 0)
    {
        time = MPI_Wtime() - time;
        printf("There are %lu prime numbers under %lu. \n", count, n);
        printf("The algorithm took %f seconds to execute \n", time);
    }

    // Finalize the MPI environment.
    MPI_Finalize();
}