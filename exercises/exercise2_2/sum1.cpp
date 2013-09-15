#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <stdio.h>
using namespace std;

#define TYPE float
TYPE sum_inv_start_low(int low, int high);
TYPE sum_inv_start_high(int low, int high);

int main(int argc, char* argv[]){
    int start=1;
    int stop=atoi(argv[1]);
    TYPE outlow=sum_inv_start_low(start,stop);
    TYPE outhigh=sum_inv_start_high(start,stop);
    printf("%15.15f, %15.15f",outlow,outhigh);
    return 0;
}


TYPE sum_inv_start_low(int low, int high){
	TYPE sum=0;
   for(int i=low;i<=high;i++){
       sum+=1.0/i;
   }
   return sum;
}
TYPE sum_inv_start_high(int low, int high){
	TYPE sum=0;
   for(int i=high;i>=low;i--){
       sum+=1.0/i;
   }
   return sum;
}
