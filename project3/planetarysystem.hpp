#ifndef PLANETARYSYSTEM_HPP
#define PLANETARYSYSTEM_HPP

//#define CALL_MEMBER_FN(object,ptrToMember)  ((object).*(ptrToMember))

#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip> 
using namespace std;

#include "celeb.hpp"
using namespace std;
	
class PlanetarySystem
{
public:
  ofstream ofile;
  double dt_2,t_2,at;
  vector<vector<double> > pos, k;
  Celeb *COM;
  int N,j;
  vector<Celeb*> nbody;
  
  ~PlanetarySystem();
  void RK4(vector<Celeb*> nbody, double dt, vector<double> t, int Nk,int n);	
  void RK4advance(Celeb *body, Celeb *com, double dt);
  
  inline double getA(Celeb *b1,Celeb *com,int axis)
	{return b1->gacc(com,axis,com->m-b1->m)*(com->m-b1->m)*(com->m-b1->m)/(com->m*com->m);}

  inline double getV(Celeb *b1,int axis){return b1->vel[axis];}
		

};

	
#endif

