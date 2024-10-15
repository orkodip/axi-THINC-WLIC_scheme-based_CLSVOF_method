//INITIALIZATION OF THE VOLUME FRACTION AND LEVEL SET FIELDS FOR A CIRCULAR INTERFACE
class INI
{
	double *X,*Y,*CX,*CY,**Fa,**Phia;	//mesh and field variables
	double Xc,Yc,rad;	//center points and radius of the circle
	double *beta_x,beta_y;	//interface thickness parameters
	double om[3],eta[3];	//normalized Gaussian weights and corresponding points
	double act_func(int i,double x,double y,double Xc,double Yc,double r);	//calculate the THINC function
	public:
		INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double xc,double yc,double *Beta_x,double Beta_y,double r);	//initialization of the input variables
		void inter();	//initial interface
		void VF();	//initial volume fractions
		void LS();	//initial level set function
};
INI::INI(double *x,double *y,double *cx,double *cy,double **fa,double **phia,double xc,double yc,double *Beta_x,double Beta_y,double r)
{
	X=x; Y=y; CX=cx; CY=cy;
	Fa=fa; Phia=phia;
	Xc=xc; Yc=yc; rad=r; beta_x=Beta_x; beta_y=Beta_y;
	om[0]=4.0/9.0; om[1]=om[2]=5.0/18.0;	//normalized Gaussian weights
	eta[0]=0.0; eta[1]=sqrt(3.0/5.0); eta[2]=-sqrt(3.0/5.0);
}
void INI::inter()
{
	const int nos=200;	//no of co-ordinates points
	double TH,pts[2][nos+1];
	for(int i=0;i<=nos;i++)	//generate circle centered at (Rc,Yc)
	{
		TH=2.0*M_PI*i/nos;
		pts[0][i]=Xc+rad*cos(TH);	//r co-ordinates
		pts[1][i]=Yc+rad*sin(TH);	//y co-ordinates
	}
	ofstream p_out("ini_inter.dat");
	p_out<<"TITLE = \"INITIAL INTERFACE\""<<endl;
	p_out<<"FILETYPE = FULL"<<endl;
	p_out<<"VARIABLES = \"r\", \"y\""<<endl;
	p_out<<"ZONE I = "<<nos+1<<", DATAPACKING = BLOCK"<<endl<<endl;
	for(int i=0;i<=nos;i++)	//print r coordinates
		p_out<<" "<<pts[0][i];
	p_out<<endl<<endl;
	for(int i=0;i<=nos;i++)	//print y coordinates
		p_out<<" "<<pts[1][i];
	p_out.close();
	cout<<"INI: INITIAL INTERFACE FILE OUTPUT SUCCESSFUL"<<endl;
}
void INI::VF()
{
	double sum,dx=X[2]-X[1];
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			sum=0.0;
			for(int l2=0;l2<3;l2++)	//Gaussian quadrature to calculate actual volume fractions
				for(int l1=0;l1<3;l1++)
					sum+=om[l2]*om[l1]*(CX[i]+0.5*dx*eta[l1])*act_func(i,eta[l1],eta[l2],(Xc-CX[i]),(Yc-CY[j]),rad);
			Fa[j][i]=0.5*(1.0+sum/CX[i]);
			if(abs(Fa[j][i])<=TRUNC) Fa[j][i]=0.0;	//truncation of insignificant values
			else if(abs(1.0-Fa[j][i])<=TRUNC) Fa[j][i]=1.0;
		}
	}
}
double INI::act_func(int i,double x,double y,double Xc,double Yc,double r)
{
	double beta,P,dx=X[2]-X[1],dy=Y[2]-Y[1];
	beta=0.5*(beta_x[i]+beta_y)/dx;
	x=0.5*dx*x; y=0.5*dy*y;
	P=r-sqrt(pow((x-Xc),2.0)+pow((y-Yc),2.0));
	return (tanh(beta*P));
}
void INI::LS()
{
	double Dc;	//distance of a point from circle center
	for(int j=1;j<=J;j++)
	{
		for(int i=1;i<=I;i++)
		{
			Dc=sqrt(pow((CX[i]-Xc),2.0)+pow((CY[j]-Yc),2.0));	//calculate distance
			Phia[j][i]=rad-Dc;	//calculate level set as the signed distance function
		}
	}
}
