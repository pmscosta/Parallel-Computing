#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#define BLOCK_LOW(id , p, n) ((id)*(n) /(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(((id) + 1), p, n) -1)
#define BLOCK_SIZE(id, p, n) ((BLOCK_LOW(((id) + 1), p, n)) - (BLOCK_LOW(id,p,n)))

void usage(void){
    std::cout << "sieve <max number>" << std::endl;
    std::cout << "<max number> range between 2 and N." <<std::endl;

}

int main (int argc, char * argv[]){
    double  elapsed_time;

    MPI_Init (&argc, &argv);

    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    int process_id;
    MPI_Comm_rank (MPI_COMM_WORLD, &process_id);

    int num_processes;
    MPI_Comm_size (MPI_COMM_WORLD, &num_processes);

    if (argc != 2){
        if(process_id == 0 ){
            usage();
            MPI_Finalize();
            exit(1);
        }
    }

    ulong   range_max= atol(argv[1]);

    ulong   sqrtn = ceil(sqrt((ulong  )range_max));

    char * pre_marked = (char *) malloc(sqrtn + 1);
    pre_marked[0] = 1;
    pre_marked[1] = 1;

    for (ulong   i = 2; i <=sqrtn; ++i){
        pre_marked[i] = 0;
    }
    ulong   pre_k=2;

    do
    {
        ulong   base= pre_k * pre_k;
        for(ulong   i =  base; i <= sqrtn; i+= pre_k)
            pre_marked[i] = 1;

        while(pre_marked[(ulong )++pre_k]);

    } while (pre_k * pre_k <= sqrtn);
    
    std::vector<int> kset;
    for(ulong   i = 3; i<=sqrtn; ++i){
        if(pre_marked[i] == 0)
            kset.push_back(i);
    }

    free(pre_marked);

    if(kset.empty()){
        std::cout << "There is 1 prime less than or equal to 2." << std::endl;
        exit(0);
    }

    ulong   low_value = 2.0 + BLOCK_LOW(process_id, num_processes, range_max - 1);
    ulong   high_value = 2.0 + BLOCK_HIGH(process_id, num_processes, range_max -1);
    ulong   block_size = BLOCK_SIZE(process_id, num_processes, range_max -1);

    if (low_value % (ulong ) 2.0 == 0){
        if(high_value % (ulong ) 2.0 == 0){
            block_size = (ulong  ) floor((ulong  ) block_size /2.0);
            high_value--;
        }
        else{
            block_size = block_size/2.0;
        }
        low_value++;
    }else{
        if(high_value % 2 == 0){
            block_size = block_size/2.0;
        }
        else{
            block_size = (ulong  )ceil((ulong  ) block_size /2.0);
        }
    }

    int temp = (range_max -1) / num_processes;

    if ((2 + temp) < (ulong  ) sqrt((ulong  ) range_max)){
        if (process_id == 0){
            std::cout << "Too many processed!" << std::endl;
            std::cout << "Process should be greater equal than sqrt(n)." << std::endl;
        }
        MPI_Finalize();
        exit(1);
    }

    char * marked = (char *) malloc (block_size);

    if(marked == NULL){
        std::cout << "Process " << process_id << " cannot allocate enough memory" << std::endl;
        MPI_Finalize();
        exit(1);
    }

    for (ulong   i =0; i < block_size;i++){
        marked[i] = 0;
    }

    ulong   first_index;
    if (process_id == 0){
        first_index=0;
    }

    ulong   kindex = 0;

    ulong   k = kset[kindex];

    ulong   count = 1;

    do{
        if(k>= low_value){
            first_index= ((k - low_value) /2) + k;
        }else if(k*k > low_value){
            first_index = (k*k - low_value) /2;
        }else{
            if (low_value % k == 0){
                first_index = 0;
            }else{
                first_index = 0;
                while((low_value + (2*first_index)) % k !=0)
                    ++first_index;
            }
        }

        for(ulong   i = first_index; i < block_size; i+=(k))
            marked[i]=1;

        k = kset[++kindex];
    }while (k*k <= range_max && kindex < (ulong  ) kset.size());

    ulong   local_count = 0;

    for (ulong   i = 0; i < block_size;i++){
        if(marked[i] == 0){
            ++local_count;
        }
    }

    free(marked); marked = 0;

    MPI_Reduce (&local_count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    elapsed_time+=MPI_Wtime();

    if (process_id == 0){
        std::cout << count << " primes found between 2 and " << range_max << std::endl;

        char st[200];
        sprintf(st, "Time: %3.3f seconds \n", elapsed_time);

        std::cout << st;
    }

    MPI_Finalize();

    return 0;
}