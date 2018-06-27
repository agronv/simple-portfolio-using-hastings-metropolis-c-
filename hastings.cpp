// Agron Velovic

#include<iostream>
#include<vector>
#include<fstream>
#include <ctime>
#include <cmath>
#include "matrix.h"
#include "DataExtraction.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Tempering function used by the Mersenne Twister.
unsigned int Temper (unsigned int N) {

   N ^= (N >> 11);
   N ^= (N << 7) & 0x9d2c5680;
   N ^= (N << 15) & 0xefc60000;
   N ^= (N >> 18);

   return N;

}

//we will use Mersienne Twister to generate random variables distributed
//uniformly between 0 and 1
double MTUniform (unsigned int seed) {
   static unsigned int X[1248], m[2], initialized = 0, k;
   unsigned int N, Y;

   // Seed the RNG when a new seed is passed or it has not yet been initialized.
   if (seed || !initialized) {
      // If no seed is specified, default is 1.
      X[0] = (seed ? seed : 1);
      // Seed X[1],X[2],...,X[623] with your favorite LCG.
      for (k = 1; k < 624; k++) {
         X[k] = 22695477 * X[k-1] + 1;
      }
      m[0] = 0; m[1] = 0x9908b0df;
      // The counter "k" is now 624.
      initialized = 1;
   }

   // Concatenate the first bit of X[k-624] with the last 31 bits of X[k-623],
   //    "|" is "bit-wise or".
   Y = (X[k-624] & 0x80000000) | (X[k-623] & 0x7fffffff);

   // Now apply the invertible linear transformation A to Y and bit-wise add X[k-227].
   X[k] = ((Y >> 1) ^ m[Y & 1] ^ X[k-227]);

   // Re-load X[0],X[1],...,X[623] as you go.
   X[k-624] = X[k];

   // Apply the tempering function.
   N = Temper (X[k]);

   // Increment the counter; shift vectors when appropriate.
   k ++;
   if (k == 1248) k = 624;

   // Now 0 <= N <= 4,294,967,295; scale it to be on the interval (0,1).
   return ( (N + 0.5) / 4294967296.0 );

}

////////////////////////////////////////////////////////////////////////////////
// This function computes elapsed computation time in seconds.//////////////////
double Time () {

   static clock_t time;
   static int initialized = 0;

   // With the first call to this function, "time" is initialized.
   if (!initialized) {
      time = clock ();
      initialized = 1;
   }

   // With subsequent calls, elapsed time since the first call is returned.
   return ((double) (clock() - time)) / CLOCKS_PER_SEC;

}

int main(){
  int size, column, i, j, TIME0;
  vector<vector<double>> returns, V, transreturns, covariance;
  vector<string> tickers;
  vector<int> w, best;

  cout<<"Enter the number of stocks you would like to use"<<endl;
  cin>>size;cout<<endl;
  column = 50;

  tickers.reserve(size);
  tickerName(tickers, size);

  returns.reserve(size);
  getVectorOfReturn(returns,tickers, size, column);

  transreturns = Transpose(returns);

  covariance = Covariance(transreturns);

  V = Scalar_Multiple(10000, covariance);

  double DeltaE, p, U, t_start, s, NewVar,OldVar, st;
  vector< vector<double> > WNew;
  int k, n;
  double E_min, T;

  s = double(size);

  // Get the temperature.
  cout<<"Enter the temperature (depends on the number of stocks given, .00005 works well for 500 stocks): "<<endl;
  cin>>T;cout<<endl;

  // Start with everything in the portfolio.  This is arbitrary. The
  //    starting portfolio should not affect the outcome. w[i] = 1 if
  //    stock "i" is in the portfolio, w[i] = 0 if not.
  for (i = 0; i < size; ++i){
     w.push_back(1);
     best.push_back(1);
  }

  vector<double> temp;
  for(i=0; i<size; ++i){
    temp.push_back(1/s);
  }
  WNew.push_back(temp);

  E_min=Var(WNew,V);

  cout<<"Minimum variance of portfolio with all stocks is "<<E_min<<endl;

  // Initialize the minimal energy and where it occurs.
  OldVar = E_min;

  // Simulate the Markov chain for specified amount of seconds.
  t_start = Time ();
  cout<<"How long would you like to run this code in seconds: "<<endl;
  cin>>TIME0; cout<<endl;
  for (n = 1; Time () - t_start < TIME0; n++) {

      // Randomly choose a "neighbor" of the current portfolio.
      // Pick i independently and uniformly from {1,...,size}.
      // a nieghbor is simply flipping one of the signs of the stocks
      // to be in our portfolio or not
      U = MTUniform (0);
      i =  int(size * U);

      if(w[i]==1){
          w[i]=0;
          s-=1.0;
      }
      else{
          w[i]=1;
          s+=1.0;
      }

      for (j=0; j<size; ++j){
          WNew[0][j]=(1.0*w[j])/s;
      }
      NewVar = Var(WNew, V);

     // Compute the change in energy associated with including/discluding
     // the ith stock
     DeltaE =  abs(NewVar-OldVar);

     // if our new state has smaller variance we will automatically move to it
     // else we would only move to our new state if p is larger than a randomly
     // generated number, the formula is given below

     if(NewVar<OldVar){
       E_min = NewVar;
       OldVar = E_min;
       for(k = 0; k<size; ++k){
         best[k] = w[k];
       }
     }
     else{
       if(T>0){
         p = exp(-DeltaE/T);
         U = MTUniform(0);
         if(p>=U){
           E_min = NewVar;
           OldVar = E_min;
           for(k = 0; k<size;++k){
             best[k] = w[k];
           }
         }
          else{
            // we will return to the originalstate if we do  not transition
            if(w[i]==1){
                w[i]=0;
                s-=1.0;
            }
            else{
                w[i]=1;
                s+=1.0;
            }
          }
           }
         }
    // we will display every step of our for loop to see how our code is progressing
     cout<<Time()<<": "<<E_min<<endl;


  } // This ends the Markov chain simulation loop.

  // we will print out the minimum variance portfolio
  cout<<endl;
  cout<<"Best found simple portfolio is: "<<endl;
  int number = 0;
  for (i = 0; i < size; ++i) {
    if (best[i]){
       cout<<tickers[i]<<": "<<1/s<<endl;
       ++number;
     }
  }
  cout<<endl;
  cout<<"Min var: "<<E_min<<endl;
  cout<<"Number of stocks in portfolio is: "<<number<<endl;
}
