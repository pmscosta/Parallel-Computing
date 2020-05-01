To install OpenMP lib run: 
`sudo apt install libomp-dev`

To install OpenMPI lib run: 
`sudo apt-get install openmpi-bin` and `sudo apt-get install libopenmpi-dev`

To compile:
* `g++ sieve.cpp -fopenmp` for the sequential and openMP
* `mpic++ sieve_mpi.cpp` or `mpic++ -O2 sieve_mpi.cpp -lm -lstdc++`

To execute:
* `./a.out <number>`
* `mpiexec ./a.out <number> <debug_flag>`      
