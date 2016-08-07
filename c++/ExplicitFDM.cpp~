#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;


double explicitCallOption(double S0,double X,double T,double r,double sigma,int iMax,int jMax)
{
 
  double S_max=2*X;
  double dS=S_max/jMax;
  double dt=T/iMax;
  // storage for the stock price and option price (old and new)
  vector<double> S(jMax+1),vOld(jMax+1),vNew(jMax+1);
 
  for(int j=0;j<=jMax;j++)
  {
    S[j] = j*dS;
  }
  //final conditions 
  for(int j=0;j<=jMax;j++)
  {
    vOld[j] = max(S[j]-X,0.);
    vNew[j] = max(S[j]-X,0.);
  }
  // loop through time 
  for(int i=iMax-1;i>=0;i--)
  {
    // apply boundary condition at S=0
    vNew[0] = 0.;
    for(int j=1;j<=jMax-1;j++)
    {
      double A,B,C;
      A=0.5*sigma*sigma*j*j*dt+0.5*r*j*dt;
      B=1.-sigma*sigma*j*j*dt;
      C=0.5*sigma*sigma*j*j*dt-0.5*r*j*dt;
      vNew[j] = 1./(1.+r*dt)*(A*vOld[j+1]+B*vOld[j]+C*vOld[j-1]);
    }
    // apply boundary condition at S=S_max
    vNew[jMax] = S[jMax] - X*exp(-r*(T-i*dt));
    // set old values to new
    vOld=vNew;
  }
  // get j* such that S_0 \in [ j*dS , (j*+1)dS ]
  int jstar;
  jstar = S0/dS;
  double sum=0.;
  // run 2 point Lagrange polynomial interpolation
  sum = sum + (S0 - S[jstar+1])/(S[jstar]-S[jstar+1])*vNew[jstar];
  sum = sum + (S0 - S[jstar])/(S[jstar+1]-S[jstar])*vNew[jstar+1];
  return sum;
}

double f(double x) {
	double pi =  4.0*atan(1.0);
	return exp(-x*x*0.5)/sqrt(2*pi);
}


// Boole's Rule (to generate normal dist)
double Boole(double StartPoint, double EndPoint, int n) {
	vector<double> X(n+1, 0.0);
	vector<double> Y(n+1, 0.0);
	double delta_x = (EndPoint - StartPoint)/double(n);
	for (int i=0; i<=n; i++) {
		X[i] = StartPoint + i*delta_x;
		Y[i] = f(X[i]);
	}
	double sum = 0;
	for (int t=0; t<=(n-1)/4; t++) {
		int ind = 4*t;
	    sum += (1/45.0)*(14*Y[ind] + 64*Y[ind+1] + 24*Y[ind+2] + 64*Y[ind+3] + 14*Y[ind+4])*delta_x;
	}
	return sum;
}

// N(0,1) cdf by Boole's Rule
double N(double x) {
	return Boole(-10.0, x, 240);
}

// Black-Scholes Call Price
double BSPrice(double S0, double K, double T, double r, double sigma, char OpType)
{
	double d = (log(S0/K) + T*(r + 0.5*sigma*sigma)) / (sigma*sqrt(T));
	double call = S0*N(d) - exp(-r*T)*K*N(d - sigma*sqrt(T));
	if (OpType=='C')
		return call;
	else
		return call - S0 + K*exp(-r*T);
}


int main()
{
 
  double S0=200.,K=150.,T=1.,r=0.05,sigma=0.3;

  int iMax=1000,jMax=100;
  cout<<"*******************"<<endl;
  cout<<"S0="<< S0<< "  K="<< K<<" Imax="<< iMax<<"  Jmax="<< jMax<< endl;
  cout<<"*******************"<<endl;
  cout <<"PDE solution = " <<explicitCallOption(S0,K,T,r,sigma,iMax,jMax) << endl;
  cout<<"*******************"<<endl;
  cout <<"Analytical solution = " << BSPrice(S0,K,T,r,sigma,'C') <<endl;
  cout<<"*******************"<<endl;

}
