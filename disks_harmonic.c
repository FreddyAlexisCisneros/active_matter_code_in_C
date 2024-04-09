#include <stdlib.h>
#include <stdio.h>
#include <time.h> 
#include <math.h>
#include <complex.h>
#include <string.h>
#include "mt19937-64.c"
#define PI (3.141592653589793)

typedef struct {
 double x,y;
 double theta;
 double new_theta;
 double Fx,Fy,E;
 double dx,dy;
} particle;
typedef struct elem elem;
struct elem
{
 int num;
 struct elem *next;
};  
typedef elem* site;
typedef struct elem_n elem_n;
struct elem_n
{
 int x;
 int y;
 struct elem_n *next;
};

typedef struct { 
 int n;
} Cells;

typedef elem_n* n_site;
site list_add(site liste, int num);
n_site list_add_n(n_site liste, int x, int y);
site list_remove(site liste, int num);
void make_neighbor_sites(n_site** N_Sites, int Lx, int Ly);
void make_neighbor_all_sites(n_site** N_Sites, int Lx, int Ly);
void get_file_name(int N, double Pe, double T, char* buffer, size_t buflen) {
 snprintf(buffer, buflen, "N_%d_Pe_%.1f_T_%.1f.txt",N,Pe,T);
}
int NEXT (int k, int L) {return (k==L-1)?0:(k+1);}
int PREV (int k, int L) {return (k==0)?(L-1):(k-1);}
double dist_PBC(double dist, double L){
  if (dist>L/2.) return dist-L;
  if (dist<-L/2.) return dist+L;
  return dist;
}
double pos_PBC(double pos, double L){
  if (pos>=L) return pos-L;
  if (pos<0) return pos+L;
  return pos;
}
int PRINT_POS(int N, int Lx, int Ly, double Dr, double Dt, int C,particle *Particles);
void InitialConditions_rd(int N, int Lx, int Ly,  site **Sites,site **Sites2, particle *Particles,particle *Particles2);
void MoveParticles(int N, int Lx, int Ly, double v, double Dr, double Dt, double kpot, double dt, site **Sites, site **Sites2, particle *Particles, particle *Particles2, double T_st, double rrate, double Switching_time, int Num_of_sims,int Gradient_region);
void add_neighbor(int k1, int k2, particle *Particles, int Lx, int Ly, double kpot, int N, double x_0, double x_1);
void add_neighbor2(int k1, int k2, particle *Particles, int Lx, int Ly, double kpot, int N, double x_0, double x_1);

int main (int argc, char *argv[]){
	
 if (argc!=14){
  printf("Input error : N Lx Ly v Dr Dt kpot dt T_st seed Switching_time Num_of_sims Gradient_region\n");
  exit(1);
 }
 
 int N;                 // The number of particles in the system.
 int Lx,Ly;             // The dimensions for the system in the x,y directions.
 int Gradient_region; 
 double v;              // The self-propulsion speed.
 double Dr,Dt,Pe;       // Rotational and translational diffusion coefficients.
 double kpot;           // The innteraction strength k.
 double dt;             // The size of the time step.
 double T_st;           // The amount of time that we will let the system take to reach a stationary-state.
 long seed;             // Seed for random number generator.
 double Switching_time; // The swithcing time.
 int Num_of_sims;       // The number of simulations that we will perform.
 clock_t start, end; 
 double cpu_time_used; 
 
 N=(int) strtol(argv[1],NULL,10);   			  // Must be an int.
 Lx=(int) strtol(argv[2],NULL,10);   			  // Must be an int.
 Ly=(int) strtol(argv[3],NULL,10); 			      // Must be an int.
 v=strtod(argv[4],NULL);            			  // Must be a double.             
 Dr=strtod(argv[5],NULL);            			  // Must be a double.
 Dt=strtod(argv[6],NULL);
 kpot=strtod(argv[7],NULL);          			  // Must be a double.
 dt=strtod(argv[8],NULL);            			  // Must be a double.
 T_st=strtod(argv[9],NULL);          			  // Must be a double.
 seed=strtol(argv[10],NULL,10);                   // Must be a double or int.
 Switching_time=strtod(argv[11],NULL);            // Must be a double.
 Num_of_sims=(int) strtol(argv[12],NULL,10);      // Must be an int.
 Gradient_region=(int) strtol(argv[13],NULL,10);
 init_genrand64(seed);   						  // Feeding the seed to the number generator.                        
 
 int Iter;                       // Will be used to iterate over in a for-loop.
 const size_t BUFLEN = 50;       // These two quantity will be used to label .txt files for outputs.
 char file_name[BUFLEN];         // These two quantity will be used to label .txt files for outputs.           
 double rrate;                   // Will be used to vary the interaction strength for the Nth particle.
 int i,j;
	 
 FILE *output;
 site **Sites;
 site **Sites2;
 
 site p;
 particle *Particles;
 particle *Particles2;
 
 int Lx_;
 Lx_ = 2*Lx + Gradient_region;
 Sites=(site **) calloc(sizeof(site *),Lx_);
 Sites2=(site **) calloc(sizeof(site *),Lx_);
 
 Particles=(particle *) malloc(sizeof(particle)*N);
 Particles2=(particle *) malloc(sizeof(particle)*N);
 
 for (i=0;i<Lx_;i++){
  Sites[i]=(site *) calloc(sizeof(site),Ly);
  Sites2[i]=(site *) calloc(sizeof(site),Ly);
 } 
 for (i=0;i<Lx_;i++) for (j=0;j<Ly;j++) Sites[i][j]=NULL;
 for (i=0;i<Lx_;i++) for (j=0;j<Ly;j++) Sites2[i][j]=NULL;

 rrate = (1./Switching_time)*dt;

 MoveParticles(N,Lx,Ly,v,Dr,Dt,kpot,dt,Sites,Sites2,Particles,Particles2,T_st,rrate,Switching_time,Num_of_sims,Gradient_region);  
 return 0;
}
void MoveParticles(int N, int Lx, int Ly, double v, double Dr, double Dt, double kpot, double dt, site **Sites,site **Sites2, particle *Particles, particle *Particles2, double T_st, double rrate, double Switching_time, int Num_of_sims, int Gradient_region){
	
 double t;
 int i,j,ii,jj,i_new,j_new;   // Site number (i,j) and new site (i_new,j_new)
 int k,k1,k2,kk;           // particle number
 particle p; 
 double x_new,y_new;    //new position
 double dx,dy;          //move
 double ampDr,ampDt;    //Amplitude of noise
 double tempp;
 int Counter;
 double Final_t;
 int Iter;
 double Pe;
 double x_0,x_1;
 int Lx_;
 double deltaV,W;
 srand((unsigned int)time(NULL)); 
 float randomNum;
 int Randy;
 
 Lx_ = Lx;
 x_0 = Lx;
 x_1 = Lx + Gradient_region;
 Lx = 2*Lx + Gradient_region;
 
 printf("N = %d Lx = %d Ly = %d v = %.3f Dr = %.3f Dt = %.3f kpot = %.3f dt = %.3f T_st = %.3f Switching_time = %.3f Num_of_sims = %d Gradient_region = %d\n",N,Lx,Ly,v,Dr,Dt,kpot,dt,T_st,Switching_time,Num_of_sims,Gradient_region); 
 
 site it1,it2;          //pointers to iterate on a site
 n_site it_site;        //to iterate on niehboring sites
 n_site **N_Sites;
 N_Sites=(n_site **) calloc(sizeof(n_site *),Lx);
 for (i=0;i<Lx;i++) {
  N_Sites[i]=(n_site *) calloc(sizeof(n_site),Ly);
 }
 
 Cells *cells;
 cells =  (Cells *) malloc(Lx*sizeof(Cells));
 for (i=0; i < Lx; i++) {
  cells[i].n = 0;
 }
 
 make_neighbor_sites(N_Sites,Lx,Ly);
 
 InitialConditions_rd(N,Lx,Ly,Sites,Sites2,Particles,Particles2);
 
 Final_t = T_st + Switching_time;
 ampDr=sqrt(2.*Dr*dt);
 ampDt=sqrt(2.*Dt*dt);
 printf("ampdt = %.3f\n",ampDt);
 double temp_pos;
 int tempyy;
 int cell;
 double Exp;
 int temp_x,temp_y;
 double fict_work;
 double next_t;
 next_t = T_st;
 double Sum;
 int mm,nn;
 int qq;
 int new_i,new_j;
 int cell_temp;
 char name[80]="";
 FILE *average_output;
 sprintf(name,"averages_for_N_%d_Lx_%d_Ly_%d_grad_%d_Dr_%.2f_Dt_%.2f_dt_%.6f.txt",N,Lx_,Ly,Gradient_region,Dr,Dt,dt);
 average_output=fopen(name,"a");
 printf("Dt = %.4f\n",Dt); 
 for (Iter=0;Iter<Num_of_sims;Iter++){
  Counter=0;
  t=0.;
  while (t < Final_t){  
	 Sum = 0.;
   t = dt*Counter;   
   Counter+=1;
   if (t > next_t){
    for (k=0;k<N;k++){
     cell=(int)Particles[k].x; 
     cells[cell].n = 0;
    }
    for (k=0;k<N;k++){
     cell=(int)Particles[k].x; 
     cells[cell].n += 1;
    }
	  for (qq=0;qq<Lx;qq++){
	   cell_temp = cells[qq].n;
	   fprintf(average_output,"%d ",cell_temp);
	  }
		fprintf(average_output,"\n");
   next_t += 2547.;
   }
   for (k=0;k<N;k++){
    Particles[k].Fx=0.;
    Particles[k].Fy=0.; 
    Particles[k].dx=0.;
    Particles[k].dy=0.; 
	  Particles[k].E=0.;
    Particles2[k].Fx=0.;
    Particles2[k].Fy=0.; 
	 	Particles2[k].E=0.; 
   }    
 
       
	 Randy = (rand() % N);
	 p=Particles[Randy];
	 i=(int)p.x;
   j=(int)p.y; 			
   for (mm=i-1;mm<i+2;mm++){ 
    for (nn=j-1;nn<j+2;nn++){ 
     temp_x = mm;
 	   temp_y = nn;
 	   if (mm < 0){
	    temp_x = Lx - 1;
	   }
	   if (nn < 0){
		  temp_y = Ly - 1;
	   }  
	   if (mm == Lx){
		  temp_x = 0;
	   }
	   if (nn == Ly){
		  temp_y = 0;
	   }	
	   it1=Sites[temp_x][temp_y];
     while(it1!=NULL){
      k1=it1->num;
      it_site=N_Sites[temp_x][temp_y];
      it2=it1->next;
      while(it2!=NULL){
       k2=it2->num;
       add_neighbor(k1,k2,Particles,Lx,Ly,kpot,N,x_0,x_1);
		   it2=it2->next;
      }
      while(it_site!=NULL){
       it2=Sites[it_site->x][it_site->y];
       while(it2!=NULL){
	      k2=it2->num;
	      add_neighbor(k1,k2,Particles,Lx,Ly,kpot,N,x_0,x_1);
	      it2=it2->next;	      
       }
       it_site=it_site->next;
      }
      it1=it1->next; 
     }
	  }
   }
	 //printf("E = %.4f\n",Particles[Randy].E);
	
	 p=Particles[Randy];
	 dx = p.Fx*dt + ampDt*gasdevMT();
   dy = p.Fy*dt + ampDt*gasdevMT();
 
 	 if (sqrt(dx*dx+dy*dy)>1.){
    printf("Error: time step is too large, dr = %lg\n",sqrt(dx*dx+dy*dy));
    exit(1);
   } 
	 Particles[Randy].dx = dx;
	 Particles[Randy].dy = dy;
   y_new=pos_PBC(p.y+dy,(double)Ly);
   // Applying fixed boundary conditions.
 	 if ( p.x + dx < 0. ||  p.x + dx > Lx ){
    x_new = p.x; 	  
   } else{ 
	  x_new = p.x + dx;
	 }
   i_new=(int)x_new;
   j_new=(int)y_new;	
	 
   Particles2[Randy].x=x_new;
   Particles2[Randy].y=y_new; 
   if (i!=i_new || j!=j_new){
    Sites2[i][j]=list_remove(Sites2[i][j],Randy);
    Sites2[i_new][j_new]=list_add(Sites2[i_new][j_new],Randy); 
   }
	 
	 if ((Particles[Randy].x < x_1) && (Particles2[Randy].x >= x_1) || (Particles[Randy].x > x_1) && (Particles2[Randy].x <= x_1)){
		//printf("x = %.4f x_new = %.4f\n",Particles[Randy].x,Particles2[Randy].x);
    for (mm=i_new-1;mm<i_new+2;mm++){ 
     for (nn=j_new-1;nn<j_new+2;nn++){
      temp_x = mm;
 	    temp_y = nn;
 	    if (mm < 0){
	     temp_x = Lx - 1;
	    }
	    if (nn < 0){
 	 	   temp_y = Ly - 1;
	    }  
	    if (mm == Lx){
		   temp_x = 0;
	    }
	    if (nn == Ly){
		   temp_y = 0;
	    }	
	    it1=Sites2[temp_x][temp_y];
      while(it1!=NULL){
       k1=it1->num;
       it_site=N_Sites[temp_x][temp_y];
       it2=it1->next;
       while(it2!=NULL){
        k2=it2->num;
        add_neighbor(k1,k2,Particles2,Lx,Ly,kpot,N,x_0,x_1);
        it2=it2->next;
       }
       while(it_site!=NULL){
        it2=Sites2[it_site->x][it_site->y];
        while(it2!=NULL){
	       k2=it2->num;
	       add_neighbor(k1,k2,Particles2,Lx,Ly,kpot,N,x_0,x_1);
	       it2=it2->next;	      
        }
        it_site=it_site->next;
       }
       it1=it1->next; 
      }
     }
    }

    deltaV = Particles2[Randy].E - Particles[Randy].E;
    W = (Particles2[Randy].Fx + Particles[Randy].Fx)*Particles[Randy].dx + (Particles2[Randy].Fy + Particles[Randy].Fy)*Particles[Randy].dy;
    Sum -= deltaV + W;
    Sum *= 0.5/Dt;
	
        if (Sum < 0.){
	   randomNum = ((float)rand())/RAND_MAX;
		 if (randomNum > exp(Sum)){
      i=(int)Particles2[Randy].x;
      j=(int)Particles2[Randy].y; 
      i_new=(int)Particles[Randy].x;
      j_new=(int)Particles[Randy].y;
      Sites2[i][j]=list_remove(Sites2[i][j],Randy);
      Sites2[i_new][j_new]=list_add(Sites2[i_new][j_new],Randy); 
      Particles2[Randy].x = Particles[Randy].x;
      Particles2[Randy].y = Particles[Randy].y; 
	   }
	  }
   }
	 
   i=(int)Particles[Randy].x;
   j=(int)Particles[Randy].y;  
   i_new=(int)Particles2[Randy].x;
   j_new=(int)Particles2[Randy].y;  
 	 Particles[Randy].x = Particles2[Randy].x;
   Particles[Randy].y = Particles2[Randy].y;
   if (i!=i_new || j!=j_new){
    Sites[i][j]=list_remove(Sites[i][j],Randy);
    Sites[i_new][j_new]=list_add(Sites[i_new][j_new],Randy); 
   }
  }
  fclose(average_output);
 }
}


void add_neighbor(int k1, int k2, particle *Particles, int Lx, int Ly, double kpot, int N, double x_0, double x_1){
 double deltax,deltay,r,amp,K1,K2;
 K1 = kpot;
 K2 = kpot;
 deltax=Particles[k2].x-Particles[k1].x;
 deltay=dist_PBC(Particles[k2].y-Particles[k1].y,(double)Ly); 
 r=sqrt(deltax*deltax+deltay*deltay);
 if ((Particles[k1].x < x_1) && (Particles[k2].x < x_1) && (r<1.)){
  if (x_0 < Particles[k1].x && Particles[k1].x < x_1){ // If the particle k1 is in the interface. 
   K1 = kpot*(x_1 - Particles[k1].x)/(x_1 - x_0);
  }
  if (x_0 < Particles[k2].x && Particles[k2].x < x_1){ // If the particle k2 is in the interface. 
   K2 = kpot*(x_1 - Particles[k2].x)/(x_1 - x_0);
  }
  amp = 0.5*(K1 + K2)*(1. - r)/r;
  Particles[k1].Fx -= amp*deltax;
  Particles[k1].Fy -= amp*deltay;
  Particles[k2].Fx += amp*deltax;
  Particles[k2].Fy += amp*deltay;
  
  Particles[k1].E += 0.5*amp*r*(1. - r);   
  Particles[k2].E += 0.5*amp*r*(1. - r);
 }
}
void InitialConditions_rd(int N, int Lx, int Ly,  site **Sites,site **Sites2, particle *Particles,particle *Particles2){
  int i,j,k;
  for (k=0;k<N;k++){
		/*
   if (k < 1127){ 
    Particles[k].x = 0.5*genrand64_real2()*(double)Lx;
    Particles[k].y = genrand64_real2()*(double)Ly;
   }
   if (k >= 1127){ 
    Particles[k].x = 0.5*Lx + 0.5*genrand64_real2()*(double)Lx;
    Particles[k].y = genrand64_real2()*(double)Ly;
   }
	 */
		/*
	 if (k == 0){
	 	Particles[k].x = 1.5;
		Particles[k].y = 1.5;
	 }
	 if (k == 1){
		// right
	 	Particles[k].x = 2.45;
		Particles[k].y = 1.5;
	 }
	 if (k == 2){
		// top
	 	Particles[k].x = 1.5;
		Particles[k].y = 2.45;
	 }
	 if (k == 3){
		// left
	 	Particles[k].x = 0.55;
		Particles[k].y = 1.5;
	 }
	 if (k == 4){
		// bottom
	 	Particles[k].x = 1.5;
		Particles[k].y = 0.55;
	 }
	 */
	 Particles[k].x = genrand64_real2()*(double)Lx;
	 Particles[k].y = genrand64_real2()*(double)Ly;
   Particles2[k].x = Particles[k].x;
   Particles2[k].y = Particles[k].y;
   i=(int)Particles[k].x;
   j=(int)Particles[k].y;
   Particles[k].theta=genrand64_real2()*2.*PI;
   Sites[i][j]=list_add(Sites[i][j],k);
   Sites2[i][j]=list_add(Sites2[i][j],k);
  }
}
int PRINT_POS(int N, int Lx, int Ly, double Dr, double Dt, int C, particle *Particles){
  int k;
  char name[50]="";
  FILE *output_pos;
  sprintf(name,"N_%d_Lx_%d_Ly_%d_Dr_%.2f_Dt_%.2f_%d.txt",N,Lx,Ly,Dr,Dt,C);
  output_pos=fopen(name,"w");
  for (k=0;k<N;k++){
	 fprintf(output_pos,"%.4f %.4f\n",Particles[k].x,Particles[k].y);
  }
  fprintf(output_pos,"\n");
  fclose(output_pos);
}
site list_add(site liste, int num){
  site new_elem=malloc(sizeof(elem));
  new_elem->num=num;
  new_elem->next=liste;
  return new_elem;
}
site list_remove(site liste, int k){
  site it,it_prec;
  it=liste;
  if (it->num==k){
    it=it->next;
    free(liste);
    return it;
  } else {
    it_prec=liste;
    it=it->next;
    while (it!=NULL){
      if (it->num==k){
	it_prec->next=it->next;
	free(it);
	return liste;
      }
      it_prec=it;
    it=it->next;
    }
  }
  printf("Erreur: element non trouve\n");
}
void make_neighbor_sites(n_site** N_Sites, int Lx, int Ly){
  int i,j;
  for (i=0;i<Lx;i++){
    for (j=0;j<Ly;j++){
      N_Sites[i][j]=list_add_n(N_Sites[i][j],NEXT(i,Lx),j); //droite
      N_Sites[i][j]=list_add_n(N_Sites[i][j],i,PREV(j,Ly)); //bas
      N_Sites[i][j]=list_add_n(N_Sites[i][j],PREV(i,Lx),PREV(j,Ly)); //bas gauche
      N_Sites[i][j]=list_add_n(N_Sites[i][j],NEXT(i,Lx),PREV(j,Ly)); //bas droite
    }
  }
}
void make_neighbor_all_sites(n_site** N_Sites_all, int Lx, int Ly){
  int i,j;
  for (i=0;i<Lx;i++){
    for (j=0;j<Ly;j++){
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],i,j);
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],NEXT(i,Lx),j);
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],i,NEXT(j,Ly));
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],i,PREV(j,Ly));
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],PREV(i,Lx),j);
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],PREV(i,Lx),PREV(j,Ly));
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],NEXT(i,Lx),NEXT(j,Ly));
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],NEXT(i,Lx),PREV(j,Ly));
      N_Sites_all[i][j]=list_add_n(N_Sites_all[i][j],PREV(i,Lx),NEXT(j,Ly));
    }
  }
}
n_site list_add_n(n_site liste, int x, int y){
  n_site new_elem=malloc(sizeof(elem_n));
  new_elem->x=x;
  new_elem->y=y;
  new_elem->next=liste;
  return new_elem;
}
