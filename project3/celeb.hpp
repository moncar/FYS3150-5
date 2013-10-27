#ifndef CELEB_HPP
#define CELEB_HPP

//#include <iostream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip> 

typedef struct {double m,x,y,vx,vy;}val;

class Celeb
{	
public:
  double pos[2],vel[2],m;//position,velocity, mass and G=
  double G;//gravitational constant using AU, yr and Mearth
  
  // declarations of various functions used by the class
  ~Celeb();
  Celeb();
  Celeb (val v);	

  double gacc(Celeb *body,int axis,double M);
  val centre(Celeb *body);
  void reset();
  void print();
  void recon(val v);
  inline double dist(Celeb *body, int axis)
  {return body->pos[axis]-pos[axis];}
  
  inline double distto(Celeb *body)//this would be sqrt(for loop of axis)
   {return sqrt(dist(body,0)*dist(body,0)+dist(body,1)*dist(body,1));}
  
  inline double dotrm(Celeb *b1, Celeb *b2,int i)//dot product m1*x1+m2*x2
   {return b1->m*b1->pos[i]+b2->m*b2->pos[i];}	
  
  inline double dotvm(Celeb *b1, Celeb *b2,int i)
   {return b1->m*b1->vel[i]+b2->m*b2->vel[i];}	

  inline double angmom(Celeb *body)//only 2d
   {return dist(body,0)*vel[1]-dist(body,1)*vel[0];}

  inline double energy(Celeb *com)//kinetic + potiential energy
   {return energykin()-G*com->m/distto(com);}//simplification assuming com.v=0
  inline double energykin()//kinetic energy
   {return 0.5*m*(vel[0]*vel[0]+vel[1]*vel[1]);}

};

#endif

