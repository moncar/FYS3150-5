#include <cstdlib> 	
#include <cstdio>	
#include <cmath>
#include <armadillo>	
#include <iostream>
#include <fstream>
#include <time.h>

using namespace std;
using namespace arma;

void gauss_elim(double a, double* b, double c, double* d, int n);
double compute_error(double* v, double*u,int n);

int main(int argc, char* argv[]){
	clock_t start, finish;

	int n = atoi(argv[1]);			//points to evaluate the problem in
	double boundary0,boundary1,ld,ud,hsq,*pd,*prhs, *pu, *pv;
	ofstream outputfile;
	///*
	
	vec d = zeros<vec>(n+2) + 2.0;	//set diagonal vector
	vec rhs = zeros<vec>(n+2);		//right-hand side vector
	vec x = linspace<vec>(0,1,n+2);	//variable vector
	
	ld = ud = -1.0;					//lower & upper diagonal values
	boundary0 = boundary1 = 0;		//boundary conditions
	hsq =  1.0/((n+1.0)*(n+1.0));	//square of step size
	
	pd = d.memptr(); prhs = rhs.memptr();	//pointers
	//Initialize rhs
	for (int i=1;i<=n;i++)
		rhs(i)=hsq*100*exp(-10.0*x(i));
	
	//set u(x)
	vec u=zeros<vec>(n+2);
	for (int i=0;i<=n+1;i++)
		u(i)= 1.0-(1.0-exp(-10.0))*x(i)-exp(-10.0*x[i]);
	
	//perform gauss elimination for tridiagonal matrix
	start=clock();
	gauss_elim(ld, pd, ud, prhs, n);
	
	//find solution; 
	vec v = rhs/d;
	cout << v(0) <<'\n'<< u(0)<< endl;
	finish = clock();
	cout<<"Tridiag: " <<(finish-start)/CLOCKS_PER_SEC<<endl;
	//*/
	/////////////////LU-decomposition comparison//////////////////
	/*
	mat A=eye<mat>(n,n);
	mat L=zeros<mat>(n,n);
	mat U=zeros<mat>(n,n);
	
	//Initialize rhs without explicit init conditions
	vec LUrhs=zeros<vec>(n);
	for (int i=0;i<=n-1;i++)
		LUrhs(i)=hsq*100*exp(-10.0*x(i));
	//Initialize the matrix A
	A(n-1,n-1)++;
	for (int i=0;i<n-1;i++){
		A(i,i)++;
		A(i,i+1)--;
		A(i+1,i)--;
	}


	start=clock();
	//Decompose A into L and U
	lu(L,U,A);
	vec vn=solve(trimatu(U),solve(trimatl(L),LUrhs));
	finish=clock();
	cout<<"LU dec: "<< (finish-start)/CLOCKS_PER_SEC<<endl;
	
	*/
	//////////computing error////////
	pu = u.memptr(); pv   = v.memptr();
	cout <<"Error: "<< compute_error(pv,pu,n) <<endl;
	
	///////////////check matrix matrix multiplication////////////////
	
	mat a=zeros<mat>(n,n);
	mat b=randu<mat>(n,n);
	mat c=randu<mat>(n,n);
	start=clock();
	for (int i=0;i<n;i++)//Row-major
		for (int k=0;k<n;k++)
			for (int j=0;j<n;j++)
				a(i,j)+=b(i,k)*c(k,j);
	finish=clock();
	cout<<"row-major: "<<(finish-start)/CLOCKS_PER_SEC<<endl;

	start=clock();			
	for (int j=0;j<n;j++)//Column major
		for (int k=0;k<n;k++)
			for (int i=0;i<n;i++)
				a(i,j)+=b(i,k)*c(k,j);
	finish=clock();
	cout<<"column-major: "<< (finish-start)/CLOCKS_PER_SEC<<endl;
	start=clock();			
	a=b*c;
	finish=clock();
	cout<<"BLAS: "<< (finish-start)/CLOCKS_PER_SEC<<endl;
	
	///////////////writing to file////////////////
	/*
	int logn= (int) log10(n);//works for our simple purposes
	ostringstream stringStream;	stringStream << "output_n"<<logn;
	string filename = stringStream.str();
	outputfile.open(filename.c_str());
	for (int i=0;i<=n+1;i++){
		outputfile<< v(i) <<' '<<u(i)<<endl;
	}
	outputfile.close();
	*/
	
	return 0;
}




void gauss_elim(double a, double* b, double c, double* r, int n){
	/* this function performs a gaussian elimination of a tridiagonal 
	 * "matrix" 
	 * where , b, the diagonal,
	 * and c,the upper tridiagonal, are n x 1 matrices (ie. vectors).
	 * b_ is the n x 1 rhs. vector (the b in Ax=b)  
	 * 
	 * n=number of points in vectors 		(int)	
	 * a: the lower tridiagonal value		(double)
	 * b: the diagonal, length=n 			(double*)
	 * c: the upper tridiagonal value		(double)
	 * d: the rhs vector(d in Ax=d), length=n
	 *
	*/
	
	//forward substitution
	double ac=a*c;
	for (int i=2;i<n+1;i++){
		b[i]-=ac/b[i-1];
		r[i]-=a*r[i-1]/b[i-1];
	}
	//backward substitution
	for (int i=n-1;i>0;i--)
		r[i]-=c*r[i+1]/b[i+1];
}	
	
double compute_error(double* v, double*u,int n){
	/*
	* compute the maximum log10 value of the relative error of
	* any index of array v of length n, compared to vector u 
	*/
	double max,errori;
	max = log10(abs((v[1]-u[1])/u[1]));
	for (int i=1; i<n+1; i++){
		errori=log10(abs((v[i]-u[i])/u[i]));
		if (errori>max) max=errori;
	}
	return max;
}
