//AXISYMMETRIC CLSVOF ALGORITHM USING WLIC-THINC SCHEME FOR VOLUME FRACTIONS AND ENO2 SCHEME FOR LEVEL SETS.
class CLSVOF:public MBASE
{
	int **tag;	//tags of interfacial cells
	double **F_act,**Phi_act;	//actual volume fraction and level set function
	double **Ft,**Phit;	//intermediate volume fraction and level set fields
	double *beta_x,beta_y;	//parameters to control slope and thickness of interface jump
	double mass_act;	//actual volume of the torus

	void updt_ghost(double **Phi);	//update the ghost cells of cell centered field
	double ENO2(int flag,int i,int j,double **Phi,double V);	//ENO2 scheme
	double UP1(int flag,int i,int j,double **Phi,double V);	//1st order upwind scheme
	double vol_frac_flux(int flag,int i,int j,double **Fa,double **Phia,double V);	//calculate the volume fraction advection flux
	void adv();	//solve the advection equations
	double LS_corr(double BETA_x,double Xc,double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02);	//enforce mass conservation to the level set function
	void reinit(int t);	//LS reinitialization algorithm
	public:
			CLSVOF(); ~CLSVOF();
			void ini(double xc,double yc,double a,double b);
			void solve(int n);	//CLSVOF advection algorithm
			double adv_err();	//calculate the advection error
			void mass_err();	//calculate mass error
			void act_vol_frac_write();	//Tecplot360 file output
			void lsvf_write(int t);	//Tecplot360 file output
			void ls_complete(int t);	//ls file output including the ghost cells
};
CLSVOF::CLSVOF()
{
	tag=new int*[J+2];
	Phi_act=new double*[J+2];
	F_act=new double*[J+1];
	Phit=new double*[J+2];
	Ft=new double*[J+2];
	beta_x=new double[I+1];
	for(int i=0;i<J+2;i++)
	{
		tag[i]=new int[I+2];
		Phi_act[i]=new double[I+2];
		Phit[i]=new double[I+2];
		Ft[i]=new double[I+2];
		if(i<J+1) F_act[i]=new double[I+1];
	}
	cout<<"CLSVOF: MEMORY ALLOCATED"<<endl;
}
CLSVOF::~CLSVOF()
{
	for(int i=0;i<J+2;i++)
	{
		delete[] tag[i];
		delete[] Phi_act[i];
		delete[] Phit[i];
		delete[] Ft[i];
		if(i<J+1) delete[] F_act[i];
	}
	delete[] tag;
	delete[] Phi_act;
	delete[] F_act;
	delete[] Phit;
	delete[] Ft;
	delete[] beta_x;
	cout<<"CLSVOF: MEMORY RELEASED"<<endl;
}
void CLSVOF::ini(double xc,double yc,double a,double b)
{
	beta_y=2.3;	//3 mesh space smoothing in z direction
	//beta_y=3.5;	//2 mesh space smoothing in z direction
	double sm=1.5;	//3 mesh space smoothing in r direction
	//double sm=1.0;	//2 mesh space smoothing in r direction
	for(int i=1;i<=I;i++)	//calculate beta in radial direction WRT CX[i]
		beta_x[i]=dx*atanh(1.0-2.0*1e-3)/(sm*(CX[i]+0.5*sm*dx));
	cout<<"CLSVOF: beta_y = "<<beta_y<<endl;

	INI ms(Xm,Ym,CX,CY,F_act,Phi_act,xc,yc,beta_x,beta_y,a);	//mTHINC algorithm is used
	ms.inter();	//export initial interface
	ms.VF();	//initial exact volume fraction is calculated here
	ms.LS();	//initial exact level set field
	mass_act=0.0;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			F[j][i]=F_act[j][i]; Phi[j][i]=Phi_act[j][i];
			mass_act+=F_act[j][i]*CX[i]*dx*dy;
		}
	}
}
void CLSVOF::updt_ghost(double **Phia)
{
	for(int j=1;j<=J;j++)	//left and right ghost nodes (Neumann bc)
	{
		Phia[j][0]=Phia[j][1];
		Phia[j][I+1]=Phia[j][I];
	}
	for(int i=1;i<=I;i++)	//bottom and top ghost nodes (Neumann bc)
	{
		Phia[0][i]=Phia[1][i];
		Phia[J+1][i]=Phia[J][i];
	}
}
double CLSVOF::ENO2(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	double plus,minus;	//plus and minus flux
	if(flag==0)	//X direction flux
	{
		plus=Phia[j][i+1]-0.5*MINMOD((Phia[j][i+1]-Phia[j][i]),(Phia[j][i+2]-Phia[j][i+1]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j][i-1]),(Phia[j][i+1]-Phia[j][i]));
	}
	else if(flag==1)	//Y direction flux
	{
		plus=Phia[j+1][i]-0.5*MINMOD((Phia[j+1][i]-Phia[j][i]),(Phia[j+2][i]-Phia[j+1][i]));
		minus=Phia[j][i]+0.5*MINMOD((Phia[j][i]-Phia[j-1][i]),(Phia[j+1][i]-Phia[j][i]));
	}
	if(V>0.0) return minus;
	else return plus;
}
double CLSVOF::UP1(int flag,int i,int j,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;
	if(flag==0)	//X direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j][i+1];
	}
	else if(flag==1)	//Y direction flux
	{
		if(V>0.0) return Phia[j][i];
		else return Phia[j+1][i];
	}
	else { cout<<"CLSVOF: ERROR IN UPWIND SCHEME!"<<endl; return 0; }
}
double CLSVOF::LS_corr(double BETA_x,double Xc,double Fa,double a_00,double a_10,double a_01,double a_11,double a_20,double a_02)
{
	double P[9],A[9];
	double beta1,om[3],eta[3];	//interface thickness parameter, Gaussian weights, and Gaussian points
	beta1=0.5*(BETA_x+beta_y)/dx;	//average of both the coefficients
	om[0]=4.0/9.0; om[1]=om[2]=5.0/18.0;	//normalized Gaussian weights
	eta[0]=0.0; eta[1]=sqrt(3.0/5.0); eta[2]=-sqrt(3.0/5.0);	//Gaussian points
	double gamma,temp,D=-1.0,C=2.0*(Fa-0.5)*Xc;	//initial guess = -1.0
	double x1,y1;
	double func,func1;	//function and its derivative
	int cnt=0;	//no of iterations of NR method
	for(int l2=0;l2<3;l2++)	//calculation of gamma
	{
		for(int l1=0;l1<3;l1++)
		{
			x1=(0.5*dx)*eta[l1]; y1=(0.5*dy)*eta[l2];
			P[l2*3+l1]=(a_00+a_10*x1+a_01*y1+a_11*x1*y1+a_20*pow(x1,2.0)+a_02*pow(y1,2.0));	//calculation of interfacial polynomial
			temp=beta1*P[l2*3+l1];
			if((l1==0)&&(l2==0)) gamma=temp;	//initialize gamma in the 1st iteration
			else if(gamma>temp) gamma=temp;
		}
	}
	gamma=-gamma+1e-6;	//final value of gamma
	for(int l2=0;l2<3;l2++)	//calculation of A_g
		for(int l1=0;l1<3;l1++)
			A[l2*3+l1]=tanh(beta1*P[l2*3+l1]+gamma);
	do	//Newton-Raphson method
	{
		func=func1=0.0;	//reinitialization
		temp=D;	//store value of previous iteration (variable is reused)
		for(int l2=0;l2<3;l2++)	//calculation of function and its derivatives
		{
			for(int l1=0;l1<3;l1++)
			{
				func+=om[l1]*om[l2]*(Xc+0.5*dx*eta[l1])*(A[l2*3+l1]+D)/(1.0+A[l2*3+l1]*D);
				func1+=om[l1]*om[l2]*(Xc+0.5*dx*eta[l1])*(1.0-pow(A[l2*3+l1],2.0))/(pow((1.0+A[l2*3+l1]*D),2.0));
			}
		}
		func-=C;
		D-=func/func1;
		cnt++;
		if(func>0.0)
		{
			cout<<"CLSVOF: ERROR IN FUNC"<<endl;
			cout<<"func = "<<func<<endl;
		}
		if(func1<0.0) cout<<"CLSVOF: ERROR IN FUNC1"<<endl;
	}
	while(abs(D-temp)>=1e-6);
	return ((atanh(D)+gamma)/beta1);
}
void CLSVOF::reinit(int t)
{
	double temp,a,b,h=dx;
	double a_00,a_01,a_10,a_11,a_20,a_02;	//interface polynomial coefficients
	double gamma=0.5*h;	//gamma is the distance parameter
//---------------------INITIALIZATION SCHEME (BASED ON ADVECTED VF FIELD)---------------------------------
	for(int j=0;j<=J+1;j++)	//reinitialize the tag values and Phit values
	{
		for(int i=0;i<=I+1;i++)
		{
			tag[j][i]=0;
			Phit[j][i]=0.0;
		}
	}
	for(int j=1;j<=J;j++)	//tagging and correction based on volume fractions
	{
		for(int i=1;i<=I;i++)
		{
			if((F[j][i]>0.2)&&(F[j][i]<0.8))	//interfacial cell
			{
				tag[j][i]=1;
				a_00=Phi[j][i];
				a_10=0.5*(Phi[j][i+1]-Phi[j][i-1])/dx;
				a_01=0.5*(Phi[j+1][i]-Phi[j-1][i])/dy;
				a_11=0.25*(Phi[j+1][i+1]-Phi[j+1][i-1]-Phi[j-1][i+1]+Phi[j-1][i-1])/(dx*dy);
				a_20=0.5*(Phi[j][i+1]-2.0*Phi[j][i]+Phi[j][i-1])/(pow(dx,2.0));
				a_02=0.5*(Phi[j+1][i]-2.0*Phi[j][i]+Phi[j-1][i])/(pow(dy,2.0));
				Phit[j][i]=Phi[j][i]+LS_corr(beta_x[i],CX[i],F[j][i],a_00,a_10,a_01,a_11,a_20,a_02);
			}
		}
	}
	for(int j=1;j<=J;j++)	//reset LS in untagged cells and correct LS values in tagged cells
	{
		for(int i=1;i<=I;i++)
		{
			if(tag[j][i]==0) Phi[j][i]=SGN(F[j][i]-0.5)*1e6;
			else if(tag[j][i]==1) Phi[j][i]=Phit[j][i];
		}
	}
	for(int j=1;j<=J;j++)	//checking in horizontal direction
	{
		for(int i=1;i<I;i++)	//excluding right boundary cells
		{
			if((SGN(Phi[j][i]*Phi[j][i+1])!=1)&&((tag[j][i]+tag[j][i+1])==0))	//LS changes sign in untagged cells
			{
				Phi[j][i]=(2.0*F[j][i]-1.0)*gamma; tag[j][i]=1;
			}
		}
	}
	for(int i=1;i<=I;i++)	//checking in vertical direction
	{
		for(int j=1;j<J;j++)	//excluding top boundary cells
		{
			if((SGN(Phi[j][i]*Phi[j+1][i])!=1)&&((tag[j][i]+tag[j+1][i])==0))	//LS changes sign in untagged cells
			{
				Phi[j][i]=(2.0*F[j][i]-1.0)*gamma; tag[j][i]=1;
			}
		}
	}
	for(int j=1;j<=J;j++)	//initialize level sets of the ghost cells(left and right boundaries)
	{
		Phi[j][0]=Phi[j][1]; tag[j][0]=tag[j][1];
		Phi[j][I+1]=Phi[j][I]; tag[j][I+1]=tag[j][I];
	}
	for(int i=0;i<=I+1;i++)	//initialize level sets of the ghost cells(bottom and top boundaries)
	{
		Phi[0][i]=Phi[1][i]; tag[0][i]=tag[1][i];
		Phi[J+1][i]=Phi[J][i]; tag[J+1][i]=tag[J][i];
	}
//---------------SOLUTION OF DISCRETE EQUATIONS(including the ghost cells)-----------------------------
	for(int sweep=1,i_ini,j_ini,di,dj;sweep<=4;sweep++)	//Gauss-Siedel sweeps
	{
		switch(sweep)	//direction of each Gauss-Siedel sweep
		{
			case 1: j_ini=0; i_ini=0;
					dj=1; di=1;
					break;
			case 2: j_ini=0; i_ini=I+1;
					dj=1; di=-1;
					break;
			case 3: j_ini=J+1; i_ini=I+1;
					dj=-1; di=-1;
					break;
			case 4: j_ini=J+1; i_ini=0;
					dj=-1; di=1;
					break;
			default: break;
		}
		for(int j=j_ini;((j>=0)&&(j<=J+1));j+=dj)	//sweep the domain in the required direction (SMART LOOPS!)
		{
			for(int i=i_ini;((i>=0)&&(i<=I+1));i+=di)
			{
				if(tag[j][i]==1) continue;	//interface cells are not updated
				if(i==0) a=Phi[j][i+1];	//left boundary
				else if(i==(I+1)) a=Phi[j][i-1];	//right boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) a=MIN2(Phi[j][i+1],Phi[j][i-1]);
					else a=MAX2(Phi[j][i+1],Phi[j][i-1]);
				}
				if(j==0) b=Phi[j+1][i];	//bottom boundary
				else if(j==(J+1)) b=Phi[j-1][i];	//top boundary
				else	//inner domain
				{
					if(SGN(Phi[j][i])==1.0) b=MIN2(Phi[j+1][i],Phi[j-1][i]);
					else b=MAX2(Phi[j+1][i],Phi[j-1][i]);
				}
				if(SGN(Phi[j][i])==1.0)	//positive viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MIN2(a,b)+h;
					else temp=0.5*(a+b+sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MIN2(Phi[j][i],temp);
				}
				else	//negative viscosity solution
				{
					if((abs(a-b)-h)>=-EPS) temp=MAX2(a,b)-h;
					else temp=0.5*(a+b-sqrt(2.0*pow(h,2.0)-pow((a-b),2.0)));
					Phi[j][i]=MAX2(Phi[j][i],temp);
				}
			}
		}
	}
}
double CLSVOF::vol_frac_flux(int flag,int i,int j,double **Fa,double **Phia,double V)
{
	if(abs(V)<=EPS) return 0.0;	//zero velocity

	int up;	//upwind cell index
	if(flag==0)	//advection in r direction
	{
		up=(V<0.0)?(i+1):i;
		if((abs(Fa[j][up])<=EPS)||(Fa[j][up]<0.0)) return 0.0;	//empty cell or over-empty cell
		else if((abs(Fa[j][up]-1.0)<=EPS)||(Fa[j][up]>1.0))	//completely filled cell or over-filled cell
			return (Xm[i]*V);
	}
	else	//advection in z direction
	{
		up=(V<0.0)?(j+1):j;
		if((abs(Fa[up][i])<=EPS)||(Fa[up][i]<0.0)) return 0.0;	//empty cell or over-empty cell
		else if((abs(Fa[up][i]-1.0)<=EPS)||(Fa[up][i]>1.0))	//completely filled cell or over-filled cell
			return V;
	}

	double Fx,Fy;	//WLIC advective flux
	double om;	//weight for WLIC scheme
	int alp,lam;	//constants of THINC scheme
	double a1,a2,Gamma2;	//interface parameters of THINC scheme
	double r_adv,Dx;	//centroid and size of advected region
	VECTOR N;	//interface normal vector

	if(flag==0)	//advection in r direction
	{
		alp=(Fa[j][up+1]>=Fa[j][up-1])?1:-1;
		a1=exp(CX[up]*beta_x[up]/dx);
		a2=exp(beta_x[up]*CX[up]*(2.0*Fa[j][up]-1.0)/(alp*dx));
		Gamma2=pow((Xm[up-1]/dx),2.0)+log((a1*a1-a1*a2)/(a1*a2-1.0))/beta_x[up];
		r_adv=CX[up]+0.5*SGN(V)*dx;	//r_face
		Dx=r_adv-sqrt(r_adv*r_adv-2.0*V*dt*r_adv);
		r_adv-=0.5*Dx;	//centroid of the advected region
		Fx=0.25*(pow((Xm[i]-Dx),2.0)-pow(Xm[i],2.0))+0.5*alp*dx*dx/beta_x[up]*log((cosh(0.5*beta_x[up]*(pow(((Xm[i]-Dx)/dx),2.0)-Gamma2)))/(cosh(0.5*beta_x[up]*(pow((Xm[i]/dx),2.0)-Gamma2))));
		Fx=-Xm[i]*V/(r_adv*Dx)*Fx;
		Fy=Xm[i]*V*Fa[j][up];
		N.x=0.5*abs(Phia[j][up+1]-Phia[j][up-1])/dx; N.y=0.5*abs(Phia[j+1][up]-Phia[j-1][up])/dy;
	}
	else	//advection in z direction
	{
		lam=(V<0.0)?0:1;
		alp=(Fa[up+1][i]>=Fa[up-1][i])?1:-1;
		a1=exp(beta_y/alp*(2.0*Fa[up][i]-1.0));
		a2=exp(beta_y);
		Gamma2=0.5/beta_y*log((a2*a2-a1*a2)/(a1*a2-1.0));	//variable is reused
		Fx=V*Fa[up][i];
		Fy=-0.5/dt*(-V*dt+alp*dy/beta_y*log((cosh(beta_y*(lam-Gamma2-V*dt/dy)))/(cosh(beta_y*(lam-Gamma2)))));
		N.x=0.5*abs(Phia[up][i+1]-Phia[up][i-1])/dx; N.y=0.5*abs(Phia[up+1][i]-Phia[up-1][i])/dy;
	}

	if((N.x+N.y)<=1e-6) om=0.5;	//WLIC scheme
	else om=N.x/(N.x+N.y);
	return (om*Fx+(1.0-om)*Fy);
}
void CLSVOF::adv()
{
	double F_adv[2],LS_adv[2];	//advection fluxes for volume fractions and LS
	int g;	//Weymouth and Yue (2010,JCP) technique
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			g=(F[j][i]>0.5)?1:0;
			if(i==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=0.0;
			}
			if(i==I)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(0,I,j,Phi,u_EW[j][I]);	//upwind scheme for right boundary flux
			}
			else	//inner domain
			{
				F_adv[1]=vol_frac_flux(0,i,j,F,Phi,u_EW[j][i]);
				LS_adv[1]=ENO2(0,i,j,Phi,u_EW[j][i]);	//ENO2 scheme for inner domain
			}
			Ft[j][i]=F[j][i]-dt/(CX[i]*dx)*(F_adv[1]-F_adv[0])+g*dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]);
			Phit[j][i]=(Phi[j][i]-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]*LS_adv[1]-Xm[i-1]*u_EW[j][i-1]*LS_adv[0]))/(1.0-dt/(CX[i]*dx)*(Xm[i]*u_EW[j][i]-Xm[i-1]*u_EW[j][i-1]));
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
		}
	}
	for(int i=1;i<=I;i++)
	{
		for(int j=1;j<=J;j++)
		{
			g=(F[j][i]>0.5)?1:0;
			if(j==1)
			{
				F_adv[0]=0.0;	//fluid does not cross the boundary
				LS_adv[0]=UP1(1,i,1,Phit,v_NS[0][i]);	//upwind scheme for bottom boundary flux
			}
			if(j==J)
			{
				F_adv[1]=0.0;	//fluid does not cross the boundary
				LS_adv[1]=UP1(1,i,J,Phit,v_NS[J][i]);	//upwind scheme for top boundary flux
			}
			else
			{
				F_adv[1]=vol_frac_flux(1,i,j,Ft,Phit,v_NS[j][i]);
				LS_adv[1]=ENO2(1,i,j,Phit,v_NS[j][i]);	//ENO2 scheme for inner domain
			}
			F[j][i]=Ft[j][i]+g*dt/dy*(v_NS[j][i]-v_NS[j-1][i])-dt/dy*(F_adv[1]-F_adv[0]);
			Phi[j][i]=Phit[j][i]*(1.0+dt/dy*(v_NS[j][i]-v_NS[j-1][i]))-dt/dy*(v_NS[j][i]*LS_adv[1]-v_NS[j-1][i]*LS_adv[0]);
			F_adv[0]=F_adv[1];	//updation for the next cell
			LS_adv[0]=LS_adv[1];
		}
	}
	updt_ghost(Phi);
}
void CLSVOF::solve(int n)
{
	adv();
	reinit(n);
}
double CLSVOF::adv_err()
{
	double num=0.0,din=0.0;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			num+=abs(F[j][i]-F_act[j][i])*CX[i]*dx*dy;
			din+=F_act[j][i]*CX[i]*dx*dy;
		}
	}
	return(num/din);
}
void CLSVOF::mass_err()
{
	double mass=0.0;
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			mass+=F[j][i]*CX[i]*dx*dy;
	cout<<"CLSVOF: RELATIVE MASS ERROR = "<<abs(mass-mass_act)/mass_act<<endl;
}
void CLSVOF::act_vol_frac_write()
{
	ofstream p_out("act_vol_frac.dat");
	p_out<<"TITLE = \"ACTUAL LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F_act\",\"Phi_act\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED)"<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F_act[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi_act[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: ACTUAL LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL"<<endl;
}
void CLSVOF::lsvf_write(int t)
{
	string fname="ls_vol_frac_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS AND VOLUME FRACTIONS\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"F\",\"Phi\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<F[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: LEVEL SETS AND VOLUME FRACTIONS FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
void CLSVOF::ls_complete(int t)
{
	string fname="ls_comlpete_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"LEVEL SETS INCLUDING GHOST CELLS\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"X\",\"Y\",\"Phi\",\"tag\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+2<<", J="<<J+2<<", DATAPACKING=BLOCK, SOLUTIONTIME="<<t*dt<<endl;
	double x=0.0-0.5*dx,y=0.0-0.5*dy;
	for(int j=0;j<=J+1;j++)	//X coordinates
	{
		for(int i=0;i<=I+1;i++)
		{
			p_out<<" "<<x;
			x+=dx;
		}
		p_out<<endl;
		x=0.0-0.5*dx;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)	//Y coordinates
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<y;
		p_out<<endl;
		y+=dy;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<Phi[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=0;j<=J+1;j++)
	{
		for(int i=0;i<=I+1;i++)
			p_out<<" "<<tag[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"CLSVOF: COMPLETE LEVEL SET FILE OUTPUT SUCCESSFUL AT n = "<<t<<endl;
}
