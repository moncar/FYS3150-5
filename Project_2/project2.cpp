#include<cstdlib>
#include<cstdio>
#include<iostream>
#include<cmath>
#include <time.h>
#include <armadillo>

using namespace std;
using namespace arma;
double maxoffdiag(mat A, int* k, int* l, int n);
void rotate(mat &A, mat &R, int k, int l, int n);
void jacobi_method(mat &A, mat &R, int n);
double offdiagsum(mat A, int n);

int main(int argc, char* argv[]){
	clock_t start, finish;
	ofstream outputfile;
	
	int n = atoi(argv[1]);			//points to evaluate the problem in
	double rhomax=atof(argv[2]);			//max value of rho
	double wr	 =atof(argv[3]);	//omega_r		
	
	mat A=zeros<mat>(n,n);
	mat R=zeros<mat>(n,n);
	double rhomin=0.0001;//needs to be nonzero for 2 electron potential
	double rhodiff=rhomax-rhomin;
	vec V=linspace<vec>(rhomin,rhomax,n);	//rho
	if (!wr)	V = V%V;		//1 electron potential=(rho^2) if wr zero
	else 		V = wr*wr*V%V+1.0/V;//2 electron interacting potential if not 
	double hsqinv=(n*n)/(rhodiff*rhodiff);	//inverted squared step size
	for (int i=0;i<n-1;i++){
		A(i,i)=2.0*hsqinv+V(i);
		A(i+1,i)=-hsqinv;
		A(i,i+1)=-hsqinv;
	}
	A(n-1,n-1)=2.0*hsqinv+V(n-1);
	
	int minint;
	if (n > 5) {
		minint = 5;}
	else {
		minint = n;}
	//Print eigenvalues using Armadillos function
	vec eigval;
	mat eigvec;
	start=clock();
	eig_sym(eigval,eigvec,A);
	finish = clock();
	cout<<"Armadillo time: " <<(finish-start)/CLOCKS_PER_SEC<<endl;
	 
	//cout<<eigvec<<endl;
	start=clock();
	jacobi_method(A,R,n);
	finish = clock();
	cout<<"Self-defined time: " <<(finish-start)/CLOCKS_PER_SEC<<endl;
	vec armamineigenval=zeros<vec>(minint); 
	//cout << A <<endl<<endl;
	for (int i=0;i<minint;i++) armamineigenval(i) = eigval(i);
	cout << "arma eigval: "<<endl<<armamineigenval(span(0,2)) << endl;
	//Print eigenvalues with self-defined functions
	vec eigenval=zeros<vec>(n);
	vec mineigenval=zeros<vec>(3);
	for (int i=0;i<n;i++) eigenval(i)=A(i,i);
	//Finding and printing the three smalles eigenvalues:
	mineigenval(0) = 99999.0;
	mineigenval(1) = 99999.0;
	mineigenval(2) = 99999.0;
	
	double s;
	for(int i=0;i<n;i++) {
		s=eigenval(i);
		if (s<mineigenval(0) and s<mineigenval(1) and s<mineigenval(2)) {
			mineigenval(2) = mineigenval(1);
			mineigenval(1) = mineigenval(0);
			mineigenval(0) = s;
		}
		else if (s<mineigenval(1) and s<mineigenval(2) and s>=mineigenval(0)) {
			mineigenval(2) = mineigenval(1);
			mineigenval(1) = s;
		}
		else if (s<mineigenval(2) and s>=mineigenval(1)) {
			mineigenval(2) = s;
		}
	}
	cout << "self,lowest eigval: "<<endl<<mineigenval << endl;
	
	cout.precision(10);
	
	cout<<"sum of off-diagonal elements: "<<offdiagsum(A,n)<<endl;
	printf("sum of off-diagonal elements: %e\n",offdiagsum(A,n));

	ostringstream stringStream;	stringStream << "output_"<<wr;
	string filename = stringStream.str();
	outputfile.open(filename.c_str());
	outputfile<< eigvec.col(0);
	outputfile.close();
	
	
}


//function which implements Jacobi's rotation algorithm LN ch 7
//int jacobi_rotation(double)
void jacobi_method(mat &A, mat &R, int n){
    //setting up eigenvector matrix
    for (int i=0;i<n;i++)   R(i,i) = 1.0;
    int k,l;
    double eps = 1e-8;
    int max_iter = (n*n*n); //double in LN
    int iterations = 0;
    double max_offdiag = maxoffdiag(A,&k,&l,n);

    while( fabs(max_offdiag) > eps && iterations < max_iter ){
        max_offdiag = maxoffdiag(A,&k,&l,n);
        rotate(A,R,k,l,n);
        iterations++;
    }
    cout << "Iterations: " << iterations << endl;
    return;
}

// Function to find max matrix element not on diagonal
double maxoffdiag(mat A, int* k, int* l, int n){
    double max=0.0;
    for (int i=0;i<n;i++){
        for (int j=i+1;j<n;j++){
            if(fabs(A(i,j)) > max){
                max = fabs(A(i,j));
                *l=i;
                *k=j;
            }
        }
    }
    return max;
}
//Function to find values of cos and sin
void rotate(mat &A, mat &R, int k, int l, int n){
    double s,c;
    if (A(k,l) != 0.0){
        double t,tau;
        //t=0;tau=0;
        tau=(A(l,l) - A(k,k))/(2*A(k,l));
        if (tau > 0)
            t=1.0/(tau+sqrt(1.0+tau*tau));//rewrite of solution of 2nd deg formula
        else
            t=-1.0/(-tau+sqrt(1.0+tau*tau));
        c=1.0/sqrt(1+t*t);
        s=c*t;

    }else{
        c=1.0;
        s=0.0;
    }
    double a_kk,a_ll,a_ik,a_il,r_ik,r_il;
    a_kk = A(k,k);//previous values of A_kk and A_ll
    a_ll = A(l,l);
    //changing matrix elements with k or l index
    A(k,k) = c*c*a_kk - 2.0*c*s*A(k,l) + s*s*a_ll;
    A(l,l) = s*s*a_kk - 2.0*c*s*A(k,l) + c*c*a_ll;
    A(k,l) = 0.0;
    A(l,k) = 0.0;

    //change remaining elements
    for (int i=0;i<n;i++){
        if (i != k && i != l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = c*a_ik - s*a_il;
            A(k,i) = A(i,k);
            A(i,l) = c*a_il + s*a_ik;
            A(l,i) = A(i,l);
        }
        //compute new eigenvectors y=Sx
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = c*r_ik - s*r_il;
        R(i,l) = c*r_il + s*r_ik;
    }
    return;
}
double offdiagsum(mat A, int n){//sum of off diagonal elements
    double sum=0.0;
    for (int i=0;i<n;i++){
        for (int j=i+1;j<n;j++){
            sum += (A(i,j));
        }
    }
    return 2*sum;//symmetric matrix, so aij+aji=2aij
}

