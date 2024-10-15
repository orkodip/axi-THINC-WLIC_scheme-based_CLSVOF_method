#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<csignal>
#define EPS 1e-12
#define TRUNC 1e-16
using namespace std;
int LC=1;	//loop condition
void control(int n)	//signal control function
{
	cout<<"PROGRAM INTERUPTED!"<<endl;
	cout<<"Enter 0 to stop: ";
	cin>>LC;
	if(LC==0) cout<<"STOP SIGNAL RECIEVED. PROGRAM WILL STOP AFTER SIMULATING CURRENT TIME-STEP."<<endl;
	else {LC=1; cout<<"SIMULATION CONTINUES."<<endl;}
}
const int I=128,J=128;	//no of grids in i and j directions
const double RADIUS=1.0,HEIGHT=1.0;	//domain size
#include "func_head.h"
#include "mbase.cpp"
#include "axi_ivf_cir.cpp"
#include "clsvof.cpp"
int main()
{
	signal(SIGINT,control);	//define signal and its function (FOR MANUALLY CONTROLLING THE PROGRAM DURING RUNTIME)
	int CNT=0;
	double dt=7.8125e-4;	//CFL = 0.1
	CLSVOF ms;
	ms.MBASE::ini(CNT,1.0,10000.0,0.07869,0.07869,1.0151,dt);	//count,rho_0,rho_1,mu_0,mu_1,sigma,dt
	ms.CLSVOF::ini(0.2,0.8,0.1,0.1);	//xc,yc,a,b
	ms.grid_write();
	ms.lsvf_write(0);
	ms.vel_ini(1);
	ms.V_ini(1);
	ms.write(0);
	for(int COUNT=CNT;(LC && (COUNT<5120));COUNT++)	//manually controlled loop
	{
		ms.CLSVOF::solve(COUNT);	//CLSVOF advection
		if((((COUNT+1)%256)==0)||(COUNT==0))
		{
			ms.lsvf_write(COUNT+1);
			ms.mass_err();
		}
	}
	ms.vel_ini(-1);
	for(int COUNT=5120;(LC && (COUNT<10240));COUNT++)	//manually controlled loop
	{
		ms.CLSVOF::solve(COUNT);	//CLSVOF advection
		if((((COUNT+1)%256)==0)||(COUNT==0))
		{
			ms.lsvf_write(COUNT+1);
			ms.mass_err();
		}
	}
	cout<<"Error due to advection = "<<ms.adv_err()<<endl;
	ms.mass_err();
	return 0;
}
