* Com optimizações, mas só funciona "bem" com inteiros (tentei mudar para unsigned long, que é o que type que devemos usar dado que é preciso ir até 2³², mas sem sucesso): https://github.com/marius92mc/sieve-of-eratosthenes-with-MPI/tree/master/src
* Sem optimizações, funciona bem com qualquer range: https://github.com/sariug/sieve-of-eratosthenes-MPI/blob/master/sieve_of_eratosthenes.cpp
* Paper: https://mmfcordeiro.files.wordpress.com/2012/10/mmfcordeiro-parallelization-of-the-sieve-of-eratosthenes.pdf 
    - o último snip de código, 1.16, parece estar bem completo, com todas as otimizações (mesmo em storage até)

* optimizations: 
    -  https://www.uio.no/studier/emner/matnat/ifi/INF3380/v12/undervisningsmateriale/inf3380-sieve.pdf
    -  https://moodle.up.pt/pluginfile.php/45745/mod_resource/content/2/CPAR_sieve.pdf