#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

//the quick 'n dirty 3-body problem:

inline double dist(double *pos){return sqrt(pos[0]*pos[0]+pos[1]*pos[1]);}
inline double dist2(double *pos1,double *pos2){return sqrt((pos1[0]-pos2[0])*(pos1[0]-pos2[0])+(pos1[1]-pos2[1])*(pos1[1]-pos2[1]));}

void derivatives(double *y, double *dydt,double *y2,double m1,double m2);
void rk4(double *pv,double h, double* vsout,double m1,double m2, double *pv2);
int main(){
	ofstream ofile;
	double Ms=333000;
	double Me=1;
	double Msa=95.2;
	double MresE=Ms+Msa;
	double MresS=Ms+Me;
	double pv[4]={1,0,0,0};//pos,vel earth
	double pvs[4]={-9.54,0,0,0};//pos,vel saturn
	pv[3]=sqrt(4*M_PI*M_PI*pv[0])/dist(pv);
	pvs[3]=-sqrt(4*M_PI*M_PI*abs(pvs[0]))/dist(pvs);
	double dt=0.001;
	int N=1e5;
	cout << " kinetic energy earth: "<<0.5*Me*(pv[2]*pv[2]+pv[3]*pv[3])<<endl;
	cout << " potential energy earth: "<<4*M_PI*M_PI/333000*Me*(-Ms-Msa/dist2(pv,pvs))<<endl;
	
	ofile.open("out");
	for (int i=0;i<N;i++){
		rk4(pv,dt,pv,Ms,Msa,pvs);
		rk4(pvs,dt,pvs,Ms,Me,pv);
		ofile << 0.0<<' '<<0.0;//to fit with plotout case 0, suns position
		for (int j=0;j<2;j++){
			ofile << setw(15) << setprecision(8) << pv[j]<<' ';
		}
		for (int j=0;j<2;j++){
			ofile << setw(15) << setprecision(8) << pvs[j]<<' ';
		}
		ofile << 0.0<<' '<<0.0;//to fit with plotout case 0, COMs position
		ofile<<endl;
	}	
	cout << " kinetic energy earth: "<<0.5*Me*(pv[2]*pv[2]+pv[3]*pv[3])<<endl;
	cout << " potential energy earth: "<<4*M_PI*M_PI/333000*Me*(-Ms-Msa/dist2(pv,pvs))<<endl;
	
	ofile.close();
	
}

void derivatives(double *y, double *dydt,double *y2,double m1,double m2)
{	//taken from program1.cpp in compphys chapter 8 
	double G=4*M_PI*M_PI/333000;	//M=solar mass = 333000 earth masses
	double r=dist(y);					//mres is here the mass of the rest of the "planets"
	double r2=dist2(y,y2);
	dydt[0]=y[2]; //derivative of x 
	dydt[1]=y[3]; //derivative of y
	dydt[2]=-G*m1*y[0]/(r*r*r)-G*m2*y[0]/(r2*r2*r2);
	dydt[3]=-G*m1*y[1]/(r*r*r)-G*m2*y[1]/(r2*r2*r2);
} // end of function derivatives  

void rk4(double *pv,double h, double* vsout,double m1, double m2,double *pv2){
	//taken from program1.cpp in compphys chapter 8 
	
	int i;
	double xh,hh,h6; 
	double *dydx, *dym, *dyt, *yt;
	//   allocate space for local vectors   
	dydx=new double [4];//k1
	dym = new double [4];//middle values k2 and k3
	dyt =  new double [4];//temporary k2, then k4
	yt =  new double [4];//temporary yt
	hh = h*0.5;
	h6 = h/6.;
	derivatives(pv,dydx,pv2,m1,m2);//k1
	
	for (i = 0; i < 4; i++) {
		yt[i] = pv[i]+hh*dydx[i];
	}
	derivatives(yt,dyt,pv2,m1,m2);	// k2   
	for (i = 0; i < 4; i++) {
		yt[i] = pv[i]+hh*dyt[i];
	}
	derivatives(yt,dym,pv2,m1,m2); //  k3 
	for (i=0; i < 4; i++) {
		yt[i] = pv[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	derivatives(yt,dyt,pv2,m1,m2);    // k4 
	//      now we upgrade y in the array yout  
	for (i = 0; i < 4; i++){
		vsout[i] = pv[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}
	delete []dym;
	delete [] dyt;
	delete [] yt;
}  
