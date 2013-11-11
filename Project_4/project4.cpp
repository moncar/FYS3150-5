#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
using namespace std;


void gauss_elim_v(double a, double* b, double c, double* r, int n);
double compute_error(double* v, double* u,int n);
int main(int argc, char* argv[]){
	double x,t,hx,ht,a,ld,ud,diag;
	int Nx,Nt,runID;
	ofstream ofile;
	if (argc<3){ cout<<"args: Nx Nt solver "<<endl;return 1;}
	runID=atoi(argv[3]);//FE=1,BE=-1,CN=0,exact=2
	Nx = atoi(argv[1]) ; Nt = atoi(argv[2]); hx=1.0/Nx; ht=0.001; a = ht/(hx*hx);
	vector<double> vn(Nx+1),vc(Nx+1),u(Nx+1),d(Nx+1);
	
	cout<<"Nx: "<<Nx<<"\thx: "<<hx<<"\talpha: "<<a<<endl;
	//exact solution:
	
	if (runID==1){//FE:
		diag=1-2*a;
		ld=ud=a;
		vn[0]=vc[0]=0;//boundary condition
		for (int i=1;i<Nx+1;i++){vc[i]=-1+hx*i;}
		cout<<"FE"<<endl;
		for (int k=0;k<Nx+1;k++);//cout<<vc[k]<<' ';
		cout<<endl;
		ofile.open("outFE");		
		for (int n=0;n<Nt;n++){
			for(int i=1;i<Nx;i++){	
				vn[i]=ld*vc[i-1]+diag*vc[i]+ud*vc[i+1];
			}
			for (int k=0;k<Nx+1;k++){
				//cout<<vn[k]<<' ';
				ofile << setw(15) << setprecision(8) << vn[k]<<' ';
			}
			//cout<<endl;
			ofile<<endl;
			for(int i=1;i<Nx;i++){	
				vc[i]=vn[i];
			}	
		}
	}
	if (runID==-1){//BE:
		diag=1.0+2*a;
		ld=ud=-a;
		vn[0]=vc[0]=0;//boundary condition
		for (int i=1;i<Nx+1;i++){vc[i]=-1+hx*i;}
		cout<<"BE"<<endl;
		for (int k=0;k<Nx+1;k++);//cout<<vc[k]<<' ';
		cout<<endl;
		ofile.open("outBE");		
		for (int i=0;i<Nx+1;i++) {d[i]=diag;}
		
		for (int n=0;n<Nt;n++){
			gauss_elim_v(ld,&d[0],ud,&vc[0],Nx);
			for (int i=0;i<Nx+1;i++){
				vn[i]=vc[i]/d[i];
				d[i]=diag;
			}
			vn[0]=0;vn[Nx]=0;
			for (int k=0;k<Nx+1;k++){
				//cout<<vn[k]<<' ';
				ofile << setw(15) << setprecision(8) << vn[k]<<' ';
			}
			//cout<<endl;
			ofile<<endl;
			
			for (int i=0;i<Nx+1;i++)
				vc[i]=vn[i];
				
		}
	}
	if (runID==0){
		diag=2.0+2.0*a;
		ld=ud=-a;
		vn[0]=vc[0]=1;
		vector<double> vtmp(Nx+1);
		for (int i=1;i<Nx+1;i++){vc[i]=-1+hx*i;}
		cout<<"CN"<<endl;
		for (int k=0;k<Nx+1;k++);//cout<<vc[k]<<' ';
		//cout<<endl;
		ofile.open("outCN");
		for (int i=0;i<Nx+1;i++) d[i]=diag;
		for (int n=0;n<Nt;n++){
			vn[Nx]=vc[Nx]=vn[0]=vc[0]=0;	
			for (int i=1;i<Nx;i++){	
				vtmp[i]=a*vc[i-1]+(2.0-2.0*a)*vc[i]+a*vc[i+1];
				if (abs(vtmp[i])<0.01 | abs(vtmp[i])>100) cout<<':'<<vtmp[i]<<'\t';
			}
			gauss_elim_v(ld,&d[0],ud,&vtmp[0],Nx);
			for (int i=0;i<Nx+1;i++){
				vn[i]=vtmp[i]/d[i];
				d[i]=diag;
				
			}
			for (int k=0;k<Nx+1;k++){
				//cout<<vn[k]<<' ';
				ofile << setw(15) << setprecision(8) << vn[k]<<' ';
			}
			//cout<<endl;
			ofile<<endl;
			for (int i=1;i<Nx;i++){
				vc[i]=vn[i];
			}
			
				
		}	
	}
	ofile.close();
	
	for (int k=0;k<Nx+1;k++)
		cout<<vc[k]<<' ';cout<<endl;
	if(runID==2) ofile.open("outEX");
	
	for (int n=0;n<Nt;n++){
		t=n*ht;
		for (int i=0;i<Nx+1;i++){ 
			x=i*hx;
			u[i]=0;
			for (int j=1;j<2000;j++){
				u[i]-=2/(j*M_PI)*sin(j*M_PI*x)*exp(-j*j*M_PI*M_PI*t);//(v(x,t)
			}
			//u[i]+=us[i];//u(x,t)=v+us
			//cout<<u[i]<<' ';
			if(runID==2) ofile << setw(15) << setprecision(8) << u[i]<<' ';
		}
		//cout<<endl;
		if(runID==2) ofile<<endl;
	}
	cout<<"exact: "<<endl;
	for (int k=0;k<Nx+1;k++)
		cout<<u[k]<<' ';cout<<endl;
	
	
	
	double err=compute_error(&vc[0],&u[0],Nx);
	cout<<"Error: "<<err<<endl;
	
}

void gauss_elim_v(double a, double* b, double c, double* r, int n){
	/* this function performs a gaussian elimination of a tridiagonal "matrix" 
	 * where , b, the diagonal,
	 * and c,the upper tridiagonal, are n+1 x 1 matrices (ie. vectors).
	 * r is the n+1 x 1 rhs. vector (the b in Ax=b)  
	 * 
	 * n=number of points in vectors 		(int)	
	 * a: the lower tridiagonal value		(double)
	 * b: the diagonal, length=n+1 			(double*)
	 * c: the upper tridiagonal value		(double)
	 * r: the rhs vector(the "b" in Ax=b), length=n+1
	 *
	*/
	
	//forward substitution
	double ac=a*c;
	for (int i=1;i<n;i++){
		b[i]-=ac/b[i-1];
		r[i]-=a*r[i-1]/b[i-1];
	}
	//backward substitution
	for (int i=n-1;i>0;i--)
		r[i]-=c*r[i+1]/b[i+1];
}
double compute_error(double* v, double* u,int n){
	/*
	* compute the maximum log10 value of the relative error of
	* any index of array v of length n, compared to vector u 
	*/
	double max,errori;
	max = log10(abs((v[1]-u[1])/u[1]));
	for (int i=0; i<n; i++){
		errori=log10(abs((v[i]-u[i])/u[i]));
		if (errori>max) max=errori;
	}
	return max;
}
