/**********************

Rails & Rong code:  Time dracts make rails of
and a two adjancent  rails make a Ladder with
planar rungs!

To clean up write a printRails routine
and findNeighbor routine

Current code is 1d Ising chain

H = - Gamma sigma^x_i - K \sigma^z_i\sigma^z_{i+1}

The Poisson process is P(t) = Gamma [ - Gamma t] 
The space percolation 1 - Exp[ -overlap@ K ]


***********************/

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <cassert>
#include <array>
#include <limits>
#include <stack>
#include <ctime>
 

  
using namespace std;

#define  Debug 0
#define  L  16



#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


struct param{
  double DeltaT;
  double Gamma;
  double K;
  int N;
  int BC;   // BC = 0 open, BC = 1 Periodic
  int site;
  int current_x;
  int current_t;
  vector<double>  time[L];
  vector<int>  spin[L];
  vector<double>  aug_time[L];
  vector<int>  aug_spin[L];
  vector<int> label[L];
};


int makeRails( param * p);  //A Initiate Configuration State
int splitRails(param *p);   //B Time Percolate
int makeCluster(param *p);  //C Connected Compoent flip spins
bool find_next_site (int * rail_x,int  * tx,  param* p); 
double makeState(param *p);  //D reconstruct  New State and return B
// put in .h file
double ran_expo(double lambda);
int  poisson(double lambda);
void printState(param p);
void printAug(param p);   // 
void printLabel(param p);

//Current Phython use spin.txt and link.txt form Make cluster. 
void  printAugSpinfile(param p);    // to spin.txt 

 

int main()
{
  param p;
  p.DeltaT = 2.0;
  p.Gamma = 1.0;  
  p.K = 10.0; 
  p.BC = 0;  // -1 COLD  same spin  +1 COLD  rails with rqndome spin
  p.site = 0;
  srand(31415); // set seed
  int Max_iter = 0;  //1024;
  double  mag = 0;
  int cluster_number =0;
  
  // Change seed to time(NULL)
#if 1
  // Change seed to time(NULL)
  int seed =  time(NULL);
  srand(seed);
  cout <<  "  time(NULL) = seed =  " << seed <<endl;
 #endif  

  cout << "Length of spin chain L =  "<<L<<"   BC = " << p.BC <<  endl;
  cout << " p.DeltaT = " <<  p.DeltaT << " p.Gamma =  " <<  p.Gamma << endl;
   
  makeRails(& p);
  cout << " After makeRails " << endl;
  if(L<33)  printState(p);
  printState(p);
  
 
#if 1
   // Test loop send data to file for Python graphics 
      splitRails(& p);
 cout << "After splitRails  time and spin data. " << endl;
  if(L<256)   printAug(p);
  
  cluster_number =  makeCluster(& p);

  cout << "After makeCluster  time and spin data.  Number of clusters =  " << cluster_number << endl;
    printAugSpinfile(p);   // spin.txt 
  cout << "Cluster Labels " << endl;
  // printAugfile(p);   // cluster.txt  not used by Python at present
  // Add cluster links links.txt
#endif

 
#if 0    
  FILE *outfile;
  outfile = fopen("mag.dat","w");
  
  for(int iter = 0; iter < Max_iter; iter++)
    { //  cout << "iter = " << iter << endl;
      splitRails(& p);
      cluster_number = makeCluster(& p);
      mag = makeState(& p);
      fprintf(outfile,"  %d      %.15f    %d   \n", iter, mag, cluster_number);
      //       cout << iter << "     " << mag << endl;
    }
  
  fclose(outfile);
#endif
  
  return 0;
  
}

 /***********
 The life time starting t0  is  P(t) = lambda exp[ - lambda (t - t0)] 
 with lambda = Gamma. Using Theorem in Documents this is periodic in
interval 0:DeltaT
  *****************/

int makeRails( param * p)
{
  double delta_t;
  int s, parity;
  double t = 0;

  
  
  for(int x = 0; x < L;x++)
    {
      //   cout <<" x = " << x << endl;
      p->time[x].clear();
      p->spin[x].clear();
      t = 0;
      s = 2*((rand()+1)%2) -1;
      parity = 1;
       cout << "x is: " << x << " and t is: " << t << endl;
      while(t <  p->DeltaT)
	{
	  // put no break in 0-the position
	  p->time[x].push_back(t);
	  p->spin[x].push_back(s);
	  parity = - parity;
	  s = -s;
	  delta_t = ran_expo(p->Gamma);
	  t = t + delta_t;
	}
      
      // If odd number start over
      if(p->time[x].size() % 2 != 1)  x--;
    }
  return  0;
}


/************************
The cases of  no break cuts  and periodic wrap around arw not implemented.

 ***********************/

int splitRails(param *p)
{
  double delta_t;
  double t;
  int s;
  double next_time = 0.0;
  
  for(int x = 0; x < L;x++)
    {
      p->aug_time[x].clear();
      p->aug_spin[x].clear();
      for(int it = 0; it < p->time[x].size(); it++)
	{
	  t =  p->time[x][it];      
	  s =  p->spin[x][it];
	  p->aug_time[x].push_back(t);
	  p->aug_spin[x].push_back(s);
	  delta_t = ran_expo(p->Gamma);
	  t = t + delta_t;
	  next_time = p->DeltaT;
	  if(it < p->time[x].size() -1)
	    next_time = p->time[x][it+1];
	  
	  while(t <   next_time  )
	    {
	     p->aug_time[x].push_back(t);
	     p->aug_spin[x].push_back(s);
	     delta_t = ran_expo(p->Gamma);
	     t = t + delta_t;
	    }
	}
    }

  return 0;
}

/***************************
Open BC is implemented only.

 ************************/

int makeCluster(param *p){
  for(int x = 0; x< L;x++)
    { p->label[x].clear();
      for(int t = 0; t< p->aug_time[x].size(); t++)
	{
	  p->label[x].push_back(0); // zero label means not been found
	}
    }
  
  int Nspins = 0 ;
  for(int x = 0; x<L; x++) Nspins +=p->label[x].size();
  
  // create the stack
  std::stack<int> stack;
  double time1x, time2x, time1y, time2y;
  int cluster_number = 0;
  double overlap = 0.0;
    int size_x, size_y;
  
  //  initialize  
  int rail_x = 0;
  int rail_y =0;
  int tx = 0;
  int ty = 0;
  p->site = -1;
  p->current_t =  0; // lexical scan of unmarked site.
  p->current_x = 0;
  double ran = -1.0;
  int spin_flip = 1;
  int max_y = -1;
  int cluster_size = 0;
  
  FILE *outfile;
  outfile = fopen("link.txt","w");
  fprintf(outfile,"# rail_x       tx    rail_y     ty       cluster     overlap      outime       intime \n");

  //emter stack with new cluster label
  while(p->site < Nspins)
    {
      cluster_number += 1;
      cluster_size = 1;
      
      cout << " cluster_number at start stack = " <<  cluster_number << endl;
      p->label[rail_x][tx] = cluster_number; 
      stack.push(rail_x); stack.push(tx);
      rand()%2 == 0 ? spin_flip = 1 : spin_flip = -1;
      p->aug_spin[rail_x][tx] = spin_flip;
      p->site++;
      p->current_x = rail_x;
      //  cout << "pushed in  x, t : (" << rail_x << " , " << tx << ") label =" << cluster_number<< endl;
      size_x =  p->aug_time[rail_x].size();
      if(p->BC && tx ==0 && size_x > 1)  p->label[rail_x][size_x-1] = cluster_number; 
      
      // find all spins in cluster           
      while(!stack.empty())
	{ 
	  tx = stack.top(); stack.pop(); rail_x = stack.top(); stack.pop();
	  time1x = p->aug_time[rail_x][tx];
	  time2x =  ((tx == size_x-1)?(p->DeltaT):( p->aug_time[rail_x][tx+1]));
	  
	  // finding neighbor in cluster ALL conditions are x-y symmetric	  
	  for( int mu = 0; mu < 2 ; mu++)
	    {
	      rail_y = (rail_x + 1- 2*mu + L)%L;
	      size_y = p->aug_time[rail_y].size();
	      max_y = (p->BC && size_y >1)? (size_y-1):(size_y);  // drops segement [t_i,p->Delta]
	      
	      for(int ty = 0; ty < max_y; ty++)
		{
		  time1y = p->aug_time[rail_y][ty];		  
		  time2y =  ((ty == size_y-1)?(p->DeltaT):( p->aug_time[rail_y][ty+1]));
		  
		  
		  /************* Begin Overlap calculation *****************/
		  
		  overlap =  MIN(time2x, time2y) - MAX(time1x,time1y);  //no ends or no cuts or t = 0  piece
		  
		  if(p->BC)
		    {
		      if((tx == 0) && (size_x > 1) && (ty != 0))  //x one side ends
			{
			  overlap += MIN(p->DeltaT, time2y) - MAX(p->aug_time[rail_x][size_x-1],time1y);
			}
		      if((ty == 0) && (size_y > 1) && (tx != 0))  //y one side ends
			{
			  overlap += MIN(p->DeltaT, time2x) - MAX(p->aug_time[rail_y][size_y-1],time1x);
			}
		      if((tx == 0) && (size_x > 1) && (ty == 0) &&(size_y > 1) ) // two sides 
			{ 
			  overlap += MIN(p->DeltaT, p->DeltaT) 
			    - MAX(p->aug_time[rail_x][size_x-1],p->aug_time[rail_y][size_y-1]);
			}
		    }
 /************* End Overlap calculation *****************/
		  
// test to add to cluster		  
		  if((overlap > 0) && (p->aug_spin[rail_x][tx] == spin_flip * p->aug_spin[rail_y][ty])
		     && (p->label[rail_y][ty] == 0))
		    {
		      ran = (double)rand()/(double)RAND_MAX;
		      if(ran < 1 - exp(- overlap*p->K))
			{
			  
			  stack.push(rail_y); stack.push(ty);// tricky when going to next rail don't know next one in time
			  cluster_size  += 1;
			  p->site++;
			  p->label[rail_y][ty] = cluster_number;
			  p->aug_spin[rail_y][ty] = spin_flip;
			  if(p->BC && ty ==0 && size_y > 1)  p->label[rail_y][size_y-1] = cluster_number; //BC on y
			  p->aug_spin[rail_y][size_y-1] = spin_flip;
			  fprintf(outfile,"  %5d      %5d    %5d     %5d     %5d    %15.10f  %15.10f  %15.10f  \n",
				  rail_x, tx, rail_y, ty,  cluster_number, overlap, MIN(time2x, time2y), MAX(time1x,time1y));
			}
		    }
		}
	    }
	  // end nearestneighbor
	}
      
      if(cluster_size ==1)  fprintf(outfile,"  %5d      %5d    %5d     %5d     %5d    %15.10f  %15.10f  %15.10f  \n",
				    rail_x, tx, rail_x, tx,  cluster_number, 0.0, MIN(time2x, time2x), MAX(time1x,time1x));
      
      // stak is empty
      bool is_found = false;
      int rail_next = -10;  int t_next = -10;
      is_found =  find_next_site(& rail_next, & t_next,  p);
      if(is_found)
	{
	  rail_x = rail_next;
	  tx = t_next;
	  p->current_x = rail_next;
	  p->current_t = t_next;
	}
      else
	{
	  fclose(outfile);
	  return cluster_number;
	}
      
    }    // start with new site.
  
  fclose(outfile);
  return 0;      // if returns 0 it has not found all spins. 
} 

// could use current as arguments in p.

  bool find_next_site(int* rail_next, int* t_next,  param* p)
  {
    int max_t, size_x;
    for(int rail = p->current_x; rail < L; rail++)
      { size_x =  p->aug_time[rail].size();
	max_t = size_x;
	for(int it =  0; it < max_t;it++)  
	  {	  
	    if( p->label[rail][it] == 0 ){
	      *rail_next = rail;
	      *t_next = it;
	      return true;
	    }
	  }
      }
    return false;
  } 
  



double  makeState(param *p){
  double mag = 0.0;
  double new_t, new_s;
  double next_t, next_s;
   for (int x = 0; x < L; x++)
    {
      p->time[x].clear();
      p->spin[x].clear();
      new_t = p->aug_time[x][0];
      new_s = p->aug_spin[x][0];
      p->time[x].push_back(new_t);
      p->spin[x].push_back(new_s);
      
      for (int t = 1; t < p->aug_time[x].size(); t++) 
	{ next_t = p->aug_time[x][t];
	  next_s = p->aug_spin[x][t];
	  if(new_s != next_s)
	  { p->time[x].push_back(new_t);
	    p->spin[x].push_back(new_s);
	    mag += next_s*(next_t - new_t);
	  }
	  new_t = next_t;
	  new_s = next_s;
	}
    }

  return mag;
}

/******************
Helper Functoins

generates a random number with a Poisson distribution. Lamda is the average number

*************/


int  poisson(double lambda)
{
  int k = 0;
  long double p = 1.0;
  double ran; 
  long double l = exp(-lambda);
  while (p >= l)
    {
      ran = (double)rand()/(double)RAND_MAX;
      p *= ran;
      k++;
    }
  return k-1;
}

/**************************
 P(t) = lambda [ - lambda t]  or lambda = Gamma

Proof:  int^t_0 dt' P(t') = 1- exp[ - lamba t])   = \int^u_0 du' = u 

==> 1 - u = exp[ - lambda t] ==> t = - log(1-u)/lambda
*********************/

double ran_expo(double lambda)
{
    double u;
    u = rand() / (RAND_MAX + 1.0);
    return -log(1- u) / lambda;
}


void printState(param p){
  cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "    "; 
      for (int t = 0; t < p.time[x].size(); t++)
	{
	  //  cout << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.time[x][t]  <<  "   " << p.spin[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;
}


void printAug(param p){

 
  cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "   "; 
      for (int t = 0; t < p.aug_time[x].size(); t++)
	{
	  //  outfile << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.aug_time[x][t]  <<  "   " << p.aug_spin[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;

}

void printStatefile(param p){

   ofstream fileOut("spin_state.txt");
   streambuf *outbuf;
   outbuf =  cout.rdbuf(fileOut.rdbuf());

   cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "   "; 
      for (int t = 0; t < p.time[x].size(); t++)
	{
	  //  outfile << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.time[x][t]  <<  "   " << p.spin[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;
  cout.rdbuf(outbuf);
  fileOut.close();

}

void printAugSpinfile(param p){

   ofstream fileOut("spin.txt");
   streambuf *outbuf;
   outbuf =  cout.rdbuf(fileOut.rdbuf());

  cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "   "; 
      for (int t = 0; t < p.aug_time[x].size(); t++)
	{
	  //  outfile << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.aug_time[x][t]  <<  "   " << p.aug_spin[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;
  cout.rdbuf(outbuf);
  fileOut.close();

}

void printLabel(param p){
  cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "   "; 
      for (int t = 0; t < p.label[x].size(); t++)
	{
	  //  cout << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.time[x][t]  <<  "   " << p.label[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;
}


void printLabelAug(param p){ 
  cout << endl;
  for (int x = 0; x < L; x++)
    { 
      cout << x << "   "; 
      for (int t = 0; t < p.label[x].size(); t++)
	{
	  //  cout << " x + L*t   "  <<  "   " <<  x + L*t  << "  " << endl ;
	  cout << p.aug_time[x][t]  <<  "   " << p.label[x][t]  <<  "   ";
	}
      
      cout << endl;
    }
  cout << endl;
}
