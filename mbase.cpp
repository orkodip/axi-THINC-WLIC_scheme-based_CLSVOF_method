/*
ABSTRACT BASE CLASS
CONVENTION USED ARE AS FOLLOWS
(r,y)->(X,Y)->(i,j)->(I,J)
V=(V_r,V_y)->(u,v)
*/
class MBASE
{
	protected:
	double *Xm,*Ym;	//the mesh variables
	double *CX,*CY;	//the node variables
	double **u,**v,**P;	//cell centered r velocity, y velocity, pressure
	double **u_EW,**v_NS;	//advection velocities at the cell faces
	double **u_EW_nm1,**v_NS_nm1;	//advection velocities at the cell faces at n-1 time level
	double **F,**Phi;	//volume fraction and level set function
	double **K,**K_EW,**K_NS;	//interface curvature at cell centers and cell faces
	double **rho,**mu;	//cell centered density and dynamic viscosity
	double rho_0,rho_1;	//density for fluid 0 and fluid 1
	double mu_0,mu_1;	//viscosity for fluid 0 and fluid 1
	double sigma;	//surface tension coefficient
	VECTOR Grav;	//acceleration due to gravity
	double dx,dy,dt;	//grid spacing, time step
	int COUNT;
	void grid_gen();	//uniform grid generator
	public:
			MBASE();	//memory allocation
			~MBASE();	//memory release
			void ini(int count,double Rho_0,double Rho_1,double Mu_0,double Mu_1,double SIG,double Dt);	//initialize the problem
			void vel_ini(int s);	//initialize velocity field for Hill vortex test
			void V_ini(int s);
			void grid_write();	//export grid for Tecplot360
			void write(int t);	//export velocity and pressure field for Tecplot360
			void write_den_vis(int t);	//export density and viscosity field for Tecplot360
			void write_extra(int n);	//export data required for calculation
			void write_bin(ofstream &p_out);	//export intermediate data in binary format
			void read_bin(ifstream &p_in);	//import intermediate data
};
MBASE::MBASE()
{
	Xm=new double[I+1];
	Ym=new double[J+1];
	CX=new double[I+1];
	CY=new double[J+1];
	P=new double*[J+2];
	u_EW=new double*[J+2];
	u_EW_nm1=new double*[J+2];
	K_EW=new double*[J+2];
	u=new double*[J+2];
	v=new double*[J+2];
	v_NS=new double*[J+1];
	v_NS_nm1=new double*[J+1];
	K_NS=new double*[J+1];
	F=new double*[J+2];
	Phi=new double*[J+2];
	rho=new double*[J+2];
	mu=new double*[J+2];
	K=new double*[J+1];
	for(int j=0;j<J+2;j++)
	{
		P[j]=new double[I+2];
		u[j]=new double[I+2];
		v[j]=new double[I+2];
		u_EW[j]=new double[I+1];
		u_EW_nm1[j]=new double[I+1];
		K_EW[j]=new double[I+1];
		F[j]=new double[I+2];
		Phi[j]=new double[I+2];
		rho[j]=new double[I+2];
		mu[j]=new double[I+2];
		if(j<J+1)
		{
			v_NS[j]=new double[I+2];
			v_NS_nm1[j]=new double[I+2];
			K_NS[j]=new double[I+2];
			K[j]=new double[I+1];
		}
	}
	cout<<"MBASE: MEMORY ALLOCATED"<<endl;
}
MBASE::~MBASE()
{
	for(int j=0;j<J+2;j++)
	{
		delete[] P[j];
		delete[] u[j];
		delete[] v[j];
		delete[] u_EW[j];
		delete[] u_EW_nm1[j];
		delete[] K_EW[j];
		delete[] F[j];
		delete[] Phi[j];
		delete[] rho[j];
		delete[] mu[j];
		if(j<J+1)
		{
			delete[] v_NS[j];
			delete[] v_NS_nm1[j];
			delete[] K_NS[j];
			delete[] K[j];
		}
	}
	delete[] Xm;
	delete[] Ym;
	delete[] CX;
	delete[] CY;
	delete[] P;
	delete[] u;
	delete[] v;
	delete[] u_EW;
	delete[] u_EW_nm1;
	delete[] K_EW;
	delete[] v_NS;
	delete[] v_NS_nm1;
	delete[] K_NS;
	delete[] F;
	delete[] Phi;
	delete[] rho;
	delete[] mu;
	delete[] K;
	cout<<"MBASE: MEMORY RELEASED"<<endl;
}
void MBASE::ini(int count,double Rho_0,double Rho_1,double Mu_0,double Mu_1,double SIG,double Dt)
{
	rho_0=Rho_0; rho_1=Rho_1;
	mu_0=Mu_0; mu_1=Mu_1;
	sigma=SIG;
	Grav.x=0.0; Grav.y=-9.8;
	dx=RADIUS/I; dy=HEIGHT/J;
	dt=Dt;
	COUNT=count;
	cout<<"MBASE: Rho_0 = "<<rho_0<<", Rho_1 = "<<rho_1<<endl;
	cout<<"MBASE: Mu_0 = "<<mu_0<<", Mu_1 = "<<mu_1<<endl;
	cout<<"MBASE: Sigma = "<<sigma<<", INTERFACE CURVATURE IS CALCULATED FROM LS"<<endl;
	cout<<"MBASE: dt = "<<dt<<endl;
	grid_gen();
}
void MBASE::grid_gen()	//uniform grid generation done here
{
	Xm[0]=Ym[0]=CX[0]=CY[0]=0.0;
	for(int i=1;i<=I;i++)	//generation of mesh and cell centers
	{
		Xm[i]=Xm[i-1]+dx;
		CX[i]=0.5*(Xm[i]+Xm[i-1]);
	}
	for(int j=1;j<=J;j++)
	{
		Ym[j]=Ym[j-1]+dy;
		CY[j]=0.5*(Ym[j]+Ym[j-1]);
	}
	cout<<"MBASE: GRID GENERATED"<<endl;
}
void MBASE::vel_ini(int s)
{
	double L=0.5;
	for(int j=1;j<=J;j++)
	{
		u_EW[j][0]=0.0;
		for(int i=1;i<=I;i++)
		{
			v_NS[0][i]=0.0;
			if(Xm[i]>0.05)
			{
				u_EW[j][i]=s*(0.05/Xm[i]+0.1*Xm[i]*(CY[j]-L)/pow(L,2.0));
				v_NS[j][i]=s*0.1*(1.0-pow(((Ym[j]-L)/L),2.0)-2.0*pow((CX[i]/L),2.0));
			}
		}
	}
}
void MBASE::V_ini(int s)
{
	double L=0.5;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			u[j][i]=s*(0.05/CX[i]+0.1*CX[i]*(CY[j]-L)/pow(L,2.0));
			v[j][i]=s*0.1*(1.0-pow(((CY[j]-L)/L),2.0)-2.0*pow((CX[i]/L),2.0));
		}
	}
}
void MBASE::grid_write()
{
	ofstream p_out("mesh.dat");
	p_out<<"TITLE = \"MESH\""<<endl;
	p_out<<"FILETYPE = GRID"<<endl;
	p_out<<"VARIABLES = \"r\",\"y\""<<endl;
	p_out<<"ZONE I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK"<<endl;
	for(int j=0;j<=J;j++)	//print r co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Xm[i];
		p_out<<endl;
	}
	p_out<<endl<<endl;
	for(int j=0;j<=J;j++)	//print y co-ordinates of mesh
	{
		for(int i=0;i<=I;i++)
			p_out<<" "<<Ym[j];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: GRID WRITE SUCCESSFULL"<<endl;
}
void MBASE::write(int t)
{
	string fname="uvp_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"FLOW AND PRESSURE FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"u\",\"v\",\"P\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<u[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<v[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<P[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_den_vis(int t)
{
	string fname="den_vis_"+to_string(t)+".dat";
	ofstream p_out(fname);
	p_out<<"TITLE = \"DENSITY AND VISCOSITY FIELD\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"rho\",\"mu\""<<endl;
	p_out<<"ZONE T=\""<<t*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1,2]=CELLCENTERED), SOLUTIONTIME="<<t*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<rho[j][i];
		p_out<<endl;
	}
	p_out<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<mu[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: DENSITY AND VISCOSITY FILE OUTPUT SUCCESSFULL AT n = "<<t<<endl;
}
void MBASE::write_extra(int n)
{
	ofstream p_out("extra.dat");
	p_out<<"TITLE = \"EXTRA DATA\""<<endl;
	p_out<<"FILETYPE = SOLUTION"<<endl;
	p_out<<"VARIABLES = \"K\""<<endl;
	p_out<<"ZONE T=\""<<n*dt<<"\", I="<<I+1<<", J="<<J+1<<", DATAPACKING=BLOCK, VARLOCATION=([1]=CELLCENTERED), SOLUTIONTIME="<<n*dt<<endl;
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
			p_out<<" "<<K[j][i];
		p_out<<endl;
	}
	p_out.close();
	cout<<"MBASE: EXTRA DATA FILE OUTPUT SUCCESSFUL"<<endl;
}
void MBASE::write_bin(ofstream &p_out)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_out.write((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_out.write((char *) &u_EW_nm1[j][i],sizeof(u_EW_nm1[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &v_NS_nm1[j][i],sizeof(v_NS_nm1[j][i]));
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			p_out.write((char *) &F[j][i],sizeof(F[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_out.write((char *) &Phi[j][i],sizeof(Phi[j][i]));
	cout<<"MBASE: INTERMEDIATE FILE OUTPUT SUCCESSFULL"<<endl;
}
void MBASE::read_bin(ifstream &p_in)
{
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &u[j][i],sizeof(u[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v[j][i],sizeof(v[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &P[j][i],sizeof(P[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_in.read((char *) &u_EW[j][i],sizeof(u_EW[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v_NS[j][i],sizeof(v_NS[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I;i++)
			p_in.read((char *) &u_EW_nm1[j][i],sizeof(u_EW_nm1[j][i]));
	for(int j=0;j<=J;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &v_NS_nm1[j][i],sizeof(v_NS_nm1[j][i]));
	for(int j=1;j<=J;j++)
		for(int i=1;i<=I;i++)
			p_in.read((char *) &F[j][i],sizeof(F[j][i]));
	for(int j=0;j<=J+1;j++)
		for(int i=0;i<=I+1;i++)
			p_in.read((char *) &Phi[j][i],sizeof(Phi[j][i]));
	cout<<"MBASE: SOLUTION INITIALIZED SUCCESSFULLY"<<endl;
}
