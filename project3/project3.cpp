#include <stdio.h>
#include <time.h>

#include <cstdlib> 	
#include <cstdio>	
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "celeb.hpp"
#include "planetarysystem.hpp"


using namespace std;

int main(int argv, char* argc[]){
  
  const int n=11;
  int Num=0; 	//number of celestial bodies in computation out of 11 possible
  vector<val> planetvalues;	
  vector<Celeb*> Nbody;
  int runID = atoi(argc[1]);
  char* filename=argc[2];
  int Nk=atof(argc[3]);//Nk-thousand points in time to be used in computation. Nk*1000*dt years computed
  double dt = 0.001;
  
  vector<double> t(Nk*1e3,0);
  for (int i=0;i<Nk*1e3;i++) t[i]=i*dt;
  //values of celestial bodies relative to earth units
  string pn[n]={"Sun","Earth","Moon","Mercury","Venus","Mars","Jupiter","Saturn","Uranus","Neptune","Pluto"};
  double m[n] ={333000, 1, 0.0123, 0.0553, 0.815, 0.107, 317.8, 95.2, 14.5, 17.1, 0.0022};// m/Mearth
  double x[n] ={0, 1, 1.00257, 0, 0, -1.52, 0, 9.58, 0, -30.05, 0};	//x/AU
  double y[n] ={0, 0, 0, 0.387, -0.723, 0, -5.20, 0, 19.20, 0, -39.24};	//y/AU
  
  /*
  double vx[n]={0, 0, 0, -1.61, 1.18, 0, 0.439, 0, -0.229, 0, 0.158}; 	//vx/(AU/yr) 
  double vy[n]={0, 1, 1-0.0344, 0, 0, -0.810, 0, 0.325, 0, -0.182, 0};	//vx/(AU/yr)
  scale the v/v_earth values by multiplying with velocity of earth in [AU]/[yr]:
  double factor= 2*M_PI;
  for (int i=0;i<11;i++){vx[i]*=factor; vy[i]*=factor;}
  */
  double vx[n]={ 1,-1,-1,-1, 1, 1, 1,-1,-1, 1, 1}; //velocity direction in x	
  double vy[n]={-1, 1, 1,-1, 1,-1, 1, 1,-1,-1, 1}; //velocity direction in y
  
  switch(runID){
  case 0:
  case 4: Num=2;break;	//Sun + Earth || Sun + Mercury
  case 1: Num=3;break;	//Sun + Earth + Moon
  case 2: Num=10;break;	//All planets + Sun
  case 3: Num=3;break;	//Sun + Earth + Saturn
  }
  planetvalues.resize(Num);
  Nbody.resize(Num);	
  int count;	
  for (int i=0,count=0;i<n & count<Num;i++){
    //skip non-relevant celestial bodies
    //if (runID==0 && (pn[i]!="Sun" && pn[i]!="Earth")) continue;
    if (runID==2 && pn[i]=="Moon") continue;
    else if (runID==3 && (pn[i]!="Sun" && pn[i]!="Earth" && pn[i]!="Saturn")) continue;
    else if (runID==4 && (pn[i]!="Sun" || pn[i]!="Mercury")) continue;
    planetvalues[count].x = x[i];
    planetvalues[count].y = y[i];
    planetvalues[count].vx=vx[i];
    planetvalues[count].vy=vy[i];
    planetvalues[count].m = m[i];
    count++;
  }
  
  for(int i=0;i<Num;i++){
    Nbody[i] = new Celeb(planetvalues[i]);
  }
  
  PlanetarySystem sys;
  sys.RK4(Nbody, dt, t, Nk, Num);	
  
  rename("out",filename);
  return 0;
}
