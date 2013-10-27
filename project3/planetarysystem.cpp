#include "planetarysystem.hpp"

void PlanetarySystem::RK4(vector<Celeb*> nbody, double dt, vector<double> t, int Nk,int n){
  /*
    Computes the position and velocity of a system of Celeb instances over a time t using 
    4th-order Runge-Kutta, via finding the centre of mass. positions are saved in a file "out"

    nbody: array of pointers to celestial bodies(Celeb instances)
    dt: step length in time
    t: array of equally spaced time values [0,1], unneccessary variable
    Nk: number of 1000 points in evaluation
    n: number of Celeb instances    
  */

  COM=new Celeb();	
  int N=1e3;//Nk=points/1000		
  
  //array over position of all values within Nk points
  pos.resize(N,vector<double>(2*(n+1),0));//N x 2n matrix with initial value 0
  dt_2=dt/2; 
  ofile.open("out");// to file with default name "out"
  
  for (int i=n-1;i>=0;i--){//celeb loop to find centre of mass
    COM->recon(COM->centre(nbody[i]));//center of mass M, X,Y,Vx,Vy
  }
  //celeb loop to adjust positions to start with COM in origin
  for (int i=0;i<n;i++){ 
    nbody[i]->pos[0]-=COM->pos[0];
    nbody[i]->pos[1]-=COM->pos[1];
  }
  //fix velocities to fit accelleration:
  for (int i=0;i<n;i++){ 
    nbody[i]->vel[0]*=sqrt(abs(nbody[i]->distto(COM)*getA(nbody[i],COM,1)));// ay=vx^2/r
    nbody[i]->vel[1]*=sqrt(abs(nbody[i]->distto(COM)*getA(nbody[i],COM,0)));// ax=vy^2/r
    nbody[i]->print();
  }
  COM->reset();



  
  for(int l=0;l<Nk;l++){//file write loop
    for (int p=0;p<N;p++){//time loop, 1000 points before every file write 		
      for (int i=n-1;i>=0;i--){//celeb loop to find centre of mass
				COM->recon(COM->centre(nbody[i]));//center of mass M, X,Y,Vx,Vy
      }//(int i=0;i<n;i++)


	//find Angmom, V, M pos Energy ,mV**2/2-gm/r
      for (int i=0;i<n;i++){ 
        cout <<nbody[i]->angmom(COM)<<'\t';
				cout <<nbody[i]->energy(COM)<<'\t';
				cout <<nbody[i]->energykin()<<'\t';
				
  		}	      
      cout<<endl;
      for (int i=0;i<n;i++){  
				PlanetarySystem::RK4advance(nbody[i],COM,dt);//x and v updated
				pos[p][2*i]=nbody[i]->pos[0];//write position to array
				pos[p][2*i+1]=nbody[i]->pos[1];
      }
      pos[p][2*n]  =COM->pos[0];
      pos[p][2*n+1]=COM->pos[1];
      //COM->print();
      //cout<<"Time: "<<t[p]<<endl;
      /*cout << COM->pos[0]<< '\t'<<COM->pos[1]<<'\t'
           <<nbody[0]->pos[0]<<'\t'<<nbody[0]->pos[1]<<"\t\t"
	   <<getA(nbody[0],COM,0)<<'\t'<<getA(nbody[0],COM,1)<<endl;
      */
      //cout << nbody[1]->distto(COM)<<'\t'<<nbody[1]->distto(nbody[0])<<endl;
      COM->reset();//reset the center of mass for each time step
    }
    //write pos to file "out"
	
    for (int i=0;i<N;i++){
      //each celeb + COM has one column for x and y 
      for(int p=0;p<2*(n+1);p++){	
				ofile << setw(15) << setprecision(8) << pos[i][p]<<' ';
      }
      
      ofile << setw(15) << setprecision(8) << t[i]<<endl;
      t[i]+=N*dt;//update for new t array
    }
  
  }
  ofile.close();
	
}
PlanetarySystem::~PlanetarySystem(){}//delete COM;}
void PlanetarySystem::RK4advance(Celeb *body, Celeb *com, double dt){//,vector<vector<double> > k
  /*
    4th order Runge-Kutta solver for one time step

    *body: array of pointers to celestial bodies(Celeb instances)
    *com:
    dt: step length in time
  */
  vector<vector<double> > k(4,vector<double>(5,0));//variable used for k in RK4	
  
  //RK4 for vx,vy,x,y
  for(int i=0;i<2;i++){j=i+2;//compute v and x simultaneously 
    k[i][1]=getA(body,com,i);//k1 dv
    k[j][1]=getV(body,i);//k1 dx
  }			

  for(int i=0;i<2;i++){j=i+2;
    k[i][0]=body->vel[i];//temporary storage of v
    k[j][0]=body->pos[i];//temporary storage of pos			
    body->vel[i]=k[i][0]+k[i][1]*dt_2;//adjusted v: v=v0+a*dt/2
    body->pos[i]=k[j][0]+k[j][1]*dt_2;//adjusted x: x=x0+v*dt/2
  }
  for(int i=0;i<2;i++){j=i+2;
    k[i][2]=getA(body,com,i);	//k2 dv
    k[j][2]=getV(body,i);	//k2 dx	
  }			
  for(int i=0;i<2;i++){j=i+2;
    body->vel[i]=k[i][0]+k[i][2]*dt_2;		
    body->pos[i]=k[j][0]+k[j][2]*dt_2;
  }
  for(int i=0;i<2;i++){j=i+2;
    k[i][3]=getA(body,com,i);//k3 dv
    k[j][3]=getV(body,i);  //k3 dx	
  }			
  for(int i=0;i<2;i++){j=i+2;
    body->vel[i]=k[i][0]+k[i][3]*dt;
    body->pos[i]=k[j][0]+k[j][3]*dt;
  }			
  for(int i=0;i<2;i++){j=i+2;
    k[i][4]=getA(body,com,i);//k4 dv
    k[i][4]=getV(body,i);//k4 dx
  }			
  for(int i=0;i<2;i++){j=i+2;			
    body->vel[i]=k[i][0]+dt/6.0*(k[i][1]+2*k[i][2]+2*k[i][3]+k[i][4]);//new velocity
    body->pos[i]=k[j][0]+dt/6.0*(k[j][1]+2*k[j][2]+2*k[j][3]+k[j][4]);//new position
  }
  //cout<<body->m<<' '<<getA(body,com,0)<<' '<<getA(body,com,1)<<endl; 
}

