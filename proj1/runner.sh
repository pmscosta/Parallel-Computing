#! /bin/bash

printf $1

if [ $1 == "C" ]; then
 g++ -O2 matrixprod.cpp -o fileout -lpapi
else
 gfortran -O2 -o fileout matrixprod.F95
fi

if [ $6 == -1 ]; then
for(( k=0; k < 3; k = k + 1 )) ;
    do
        for (( i="$3"; i <= "$4"; i = i + "$5")) ;
            do  
                for(( j=128; j <= 512; j = j * 2)) ;
                    do ./fileout "$2" "$i" "$i" "$j";
                done;
        done;
    done;
else
for (( i="$3"; i <= "$4"; i = i + "$5")) ;
    do ./fileout "$2" "$i" "$i" "$6";
done;
fi