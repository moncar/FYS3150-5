/*
Float to binary 



use modula operator on float n
output is lowest digit in 
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
using namespace std;
int main(int argc, char* argv[]){
	double nf=atof(argv[1]);
	int ni=(int) nf;
	double dec=nf-ni;
	
	const int MAX_NUM=20;
	int counter=0;
	int* p = &counter;
	int array[2*MAX_NUM];
	for(int i=0;i<2*MAX_NUM;i++){array[i]=0;}
	
	if(ni>pow(2,MAX_NUM)-1){
		cout << "too big number"<<endl;	
	}
	
	while(ni>0 && counter<MAX_NUM){
		array[MAX_NUM-counter-1]=ni%2;
		(*p)++;
		ni/=2;
	}
	
	*p=10;
	while(dec!=0 && counter<2*MAX_NUM){
		dec*=2;
		array[counter]=  dec >= 1 ? 1 : 0;
		dec-= dec >= 1 ? 1 : 0;	
		(*p)++;	 		
	}
	for(int i=0;i<MAX_NUM;i++)
		cout << array[i];
	cout <<',';
	for(int i=MAX_NUM;i<2*MAX_NUM;i++)
		cout << array[i];
	 cout <<endl;
	return 0;	
} 
