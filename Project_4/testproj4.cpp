#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;


void gauss_elim_v(double a, double* b, double c, double* r, int n);
void gauss_elim_v2(double a, double* b, double c, double* r, int n);
int main(int argc, char* argv[]){
	double x,t,hx,ht,a,ld,ud,diag,lambdasq;
	int Nx,Nt,runID;
	ofstream ofile;
	if (argc<3){ cout<<"args: Nx Nt solver "<<endl;return 1;}
	runID=atoi(argv[3]);//FE=1,BE=-1,CN=0,exact=2
	Nx = atoi(argv[1]) ; Nt = atoi(argv[2]); hx=1.0/Nx; ht=0.01; a = ht/(hx*hx);
	//pc = vc.memptr(); pn = vn.memptr();	//pointers
	//gauss_elim(ld,pc,ud,pn,Nx-1);
	//d=zeros<vec>(Nx); 
	vector<double> us(Nx+1),vn(Nx+1),vc(Nx+1),u(Nx+1),d(Nx+1);
	
	cout<<"Nx: "<<Nx<<"\thx: "<<hx<<"\talpha: "<<a<<endl;
	for (int i=0;i<Nx+1;i++) us[i]+=1.0-i*hx;
	for (int n=0;n<Nt;n++){
		t=n*ht;
		for (int i=0;i<Nx+1;i++){ 
			x=i*hx;
			u[i]=0;
			for (int j=1;j<2000;j++){
				u[i]-=2/(j*M_PI)*sin(j*M_PI*x)*exp(-j*j*M_PI*M_PI*t);
			}
			u[i]+=us[i];
			cout<<u[i]<<' ';
		}
		cout<<endl;
		
	}
	if (runID==1){//FE:
		diag=1-2*a;
		ld=ud=a;
		vn[0]=vc[0]=1;//boundary condition
		cout<<"FE"<<endl;
		for (int k=0;k<Nx+1;k++)
			cout<<vn[k]<<' ';
		cout<<endl;
		ofile.open("outFE");		
		for (int n=0;n<Nt;n++){
			for(int i=1;i<Nx;i++){	
				vn[i]=ld*vc[i-1]+diag*vc[i]+ud*vc[i+1];
			}
			for (int k=0;k<Nx+1;k++)
				cout<<vn[k]<<' ';
			cout<<endl;
			for(int i=1;i<Nx;i++){	
				vc[i]=vn[i];
			}	
		}
	}
	if (runID==-1){//BE:
		diag=1+2*a;
		ld=ud=-a;
		vn[0]=vc[0]=1;//boundary condition
		cout<<"BE"<<endl;
		for (int k=0;k<Nx+1;k++)
			cout<<vn[k]<<' ';
		cout<<endl;
		ofile.open("outBE");		
		for (int i=0;i<Nx+1;i++) {d[i]=diag;}
		
		for (int n=0;n<Nt;n++){
			gauss_elim_v(ld,&d[0],ud,&vc[0],Nx+1);
			for (int i=0;i<Nx+1;i++){
				vn[i]=vc[i]/d[i];
				d[i]=diag;
			}
			vn[0]=1;vn[Nx]=0;
			for (int k=0;k<Nx+1;k++)
				cout<<vn[k]<<' ';
			cout<<endl;
			
			for (int i=0;i<Nx+1;i++)
				vc[i]=vn[i];
				
		}
	}
	if (runID==0){
		diag=2+2*a;
		ld=ud=-a;
		vn[0]=vc[0]=1;
		
		cout<<"CN"<<endl;
		for (int k=0;k<Nx+1;k++)
			cout<<vn[k]<<' ';
		cout<<endl;
		ofile.open("outCN");
		for (int i=0;i<Nx+1;i++) d[i]=diag;
		for (int n=0;n<Nt;n++){
			vn[Nx]=vc[Nx]=0;	
			for (int i=1;i<Nx;i++){	
				vn[i]=a*vc[i-1]+(2-2*a)*vc[i]+a*vc[i+1];//tmp use
			}
			gauss_elim_v(ld,&d[0],ud,&vn[0],Nx-1);
			for (int i=1;i<Nx;i++){
				vn[i]=vn[i]/d[i];
				d[i]=diag;
			}
			for (int k=0;k<Nx+1;k++)
				cout<<vn[k]<<' ';
			cout<<endl;
			for (int i=1;i<Nx-1;i++){
				vc[i]=vn[i];
			}
			
				
		}	
	}
	ofile.close();	
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
	void gauss_elim_v2(double a, double* b, double c, double* r, int n){
	/* this function performs a gaussian elimination of a tridiagonal "matrix" 
	 * where , b, the diagonal,
	 * and c,the upper tridiagonal, are n x 1 matrices (ie. vectors).
	 * r is the n x 1 rhs. vector (the b in Ax=b)  
	 * 
	 * n=number of points in vectors 		(int)	
	 * a: the lower tridiagonal value		(double)
	 * b: the diagonal, length=n 			(double*)
	 * c: the upper tridiagonal value		(double)
	 * r: the rhs vector(the "b" in Ax=b), length=n
	 *
	*/
	//forward substitution
	double ac=a*c;
	for (int i=1;i<n;i++){//b u r
		b[i]-=ac/b[i-1];
		r[i]-=a*r[i-1]/b[i-1];
	}
	//backward substitution
	for (int i=n-2;i>0;i--)
		r[i]-=c*r[i+1]/b[i+1];
}	
