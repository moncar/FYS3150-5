#include "celeb.hpp"
//celestial body class

//constructor	
Celeb::Celeb(val v){
  G=4.0*M_PI/333000.0;
  pos[0]=v.x;  pos[1]=v.y; 
  vel[0]=v.vx; vel[1]=v.vy;
  m=v.m;
}
Celeb::Celeb(){//alternative constructor
  pos[0]=pos[1]=vel[0]=vel[1]=m=0;
}


Celeb::~Celeb() {}// destructor

void Celeb::reset(){pos[0]=0; pos[1]=0; vel[0]=0; vel[1]=0; m=0;}//reset values of celeb instance to 0
void Celeb::recon(val v){pos[0]=v.x; pos[1]=v.y; vel[0]=v.vx; vel[1]=v.vy; m=v.m;}//set values equal to val v

void Celeb::print(){std::cout
	<<" x: "  << std::setw(15) << std::setprecision(8)<<pos[0]
	<<" y: "  << std::setw(15) << std::setprecision(8)<<pos[1]
	<<" vx: " << std::setw(15) << std::setprecision(8)<<vel[0]
	<<" vy: " << std::setw(15) << std::setprecision(8)<<vel[1]
	<<" m: "  << std::setw(15) << std::setprecision(8)<<m<<std::endl;}

double Celeb::gacc(Celeb *body,int axis,double M){//gravitational acceleration a=-G*M*x/r^3//
  double r = distto(body);
  return G*M*(dist(body,axis))/(r*r*r);
}
	
val Celeb::centre(Celeb *body){//Center of mass system values
  //breaks down if both cbodies have zero mass
  val com; 
  com.m=body->m+m;
  com.x =dotrm(body,this,0)/com.m;
  com.y =dotrm(body,this,1)/com.m;
  com.vx=dotvm(body,this,0)/com.m;
  com.vy=dotvm(body,this,1)/com.m;
  return com;
}
