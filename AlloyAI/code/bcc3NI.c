#include <stdio.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#define KB 8.61173305e-5
#define ITERATIONMAX 1500
#define ITERATIONMIN 50
#define NoPROBCONDITION 5
/*=====================================phase_type==============================================
						phase_type       			phase	
==============================================================================================
							1						A2
							2						B2
							3						D03
							5						B32
							6						F43m
							7						L21
							8                       L21'
 ================================================================================================*/
/*=====================================phase_type==============================================
						phase_type       			phase	
==============================================================================================
							1						A2
							2						B2
							3						D03(I)(i=j)
							4						D03(II)(k=l)
							5						B32
							6						F43m
							7						A2, all 0.5
							9 						L21	
							
							3,4,9 					D03	
							7 2 	                B2
 ================================================================================================*/
void bccNI3ele(double *pairen, double *inipro, double tem, double *elemu, double Tau, double *returnvalue,double *p4bx);
void bcc3multiproNI( double *paireE,int phase,double tem,double mu1,double mu2,double Tau, double *returnvalue,double *p4bx);
void bcc3phasetypedata(int s,int region, char data[10] ,char boundary[16],double *paireE, double tem, double *mu1,double *mu2,double Tau,double Granddiff,double Mudiff);
int bcc3judgetype(double *p);

//the pairen parameter input order is w1AB, w2AB, wAbAB, wABBB, w1AC,w2AC,.. W1BC, w2BC,, ...
//   ABC is FeCrAl
int main()     
{
	double paireE[12]={0.048040068,-0.023058161,0.010086594,0.009186436,-0.09221762,-0.038328291,0.009263967,0.026809306,-0.018562867,0.00035774,-0.013367108,0.018484173};
	/*========================1200K=====================================
	double tem=1200;
	
	double mu1[3]={-0.92,0.55,0}; 
	mu1[2]=(mu1[1]-mu1[0])/183;
    double mu2[3]={-1.3,1.64,0.6};     
	mu2[2]=(mu2[1]+1.05)/183;
	mu2[0]=-1.1;
	mu2[1]=-0.92;

    //phasetypedata(0,0,"bcc33.dat","bccbound33.dat",paireE,tem,mu1,mu2,1e-8,1e-5,1e-5);
	phasetypedata(1,0,"bcc33.dat","bccbound33.dat",paireE,tem,mu1,mu2,1e-8,1e-5,1e-5);
	*/
	
	
	
	double tem=2000;
	double mu1[3]={-1.3,0.8,0}; 
	mu1[2]=(mu1[1]-mu1[0])/81;
    double mu2[3]={-1.2,1.8,0.6};     
	mu2[2]=(mu2[1]-mu2[0])/81;
    bcc3phasetypedata(0,0,"bccd.dat","bccf3z.dat",paireE,tem,mu1,mu2,1e-8,1e-5,1e-5);
    
    /*
    double tem1=1400;
	double mu11[3]={-1,0.5,0}; 
	mu11[2]=(mu11[1]-mu11[0])/81;
    double mu12[3]={-1.1,1.65,0.6};     
	mu12[2]=(mu12[1]-mu12[0])/81;
    bcc3phasetypedata(0,0,"bccd.dat","bccf32.dat",paireE,tem1,mu11,mu12,1e-8,1e-5,1e-5);
	*/
	
	//sound
	Beep(784, 700);
	Beep(880, 800);
	Beep(932, 900);
	Beep(1046, 1000);
	Beep(1175, 1100);
	
	return 0;
}

void bccNI3ele(double *pairen, double *inipro, double tem, double *elemu, double Tau, double *returnvalue,double *p4bx)
{	
    //AT GIVEN temperature and chemical potential,
	//give the guess pairprobablilities to Iterate final probabilites
	//returnvalue array, tem,mu,8 point pro, NIF, energy;
	
	//paireE 1d 
	//*inipro initial point probabilities
	
	/*declaration of variables*/
	double beta=1.0/(tem*KB);
	double tau=Tau;
	double epsilon;
	double mu,expe,energyx;
	
	int i,j,k,l;
	double paireE1[3][3]={0};
	double paireE2[3][3]={0};
	double manybodyE[3][3][3][3]={0};
	double pairmu[3]={0};
//	
//	FILE *fp;
//	fp=fopen("bccenergy.dat","w");
	/*declaration of NI variables*/
	int F=0; 	//the number of iteratin, started from 0. 
	double p3aijk[3][3][3],p3bijk[3][3][3];
	double p3aijl[3][3][3],p3bijl[3][3][3];
	double p3aikl[3][3][3],p3bikl[3][3][3];
	double p3ajkl[3][3][3],p3bjkl[3][3][3];
	double p2aij[3][3],p2bij[3][3];
	double p2aik[3][3],p2bik[3][3];
	double p2ail[3][3],p2bil[3][3];
	double p2ajk[3][3],p2bjk[3][3];
	double p2ajl[3][3],p2bjl[3][3];
	double p2akl[3][3],p2bkl[3][3];
	double p1ai[3],p1bi[3];
	double p1aj[3],p1bj[3];
	double p1ak[3],p1bk[3];
	double p1al[3],p1bl[3];
	
	double p1a[3];
	double p1b[3]={0};
	double p2a[3][3];
	double p2b[3][3]={0};
	double p3a[3][3][3];
	double p3b[3][3][3]={0};
	double p4b[3][3][3][3];
    double p4a[3][3][3][3]={0};
	double deltaP[117];
	double delta=1;
	
	//initial p1a,p2a,p3a.
	for(i=0;i<3;i++) //has been changed 
	{
		p1ai[i]=inipro[i];
		p1aj[i]=inipro[3+i];
		p1ak[i]=inipro[6+i];
		p1al[i]=inipro[9+i];
	}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			p2aij[i][j]=p1ai[i]*p1aj[j];
			p2aik[i][j]=p1ai[i]*p1ak[j];
			p2ail[i][j]=p1ai[i]*p1al[j];
			p2ajk[i][j]=p1aj[i]*p1ak[j];
			p2ajl[i][j]=p1aj[i]*p1al[j];
			p2akl[i][j]=p1ak[i]*p1al[j];
	
		}
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
			{
				p3aijk[i][j][k]=p1ai[i]*p1aj[j]*p1ak[k];
				p3aijl[i][j][k]=p1ai[i]*p1aj[j]*p1al[k];
				p3aikl[i][j][k]=p1ai[i]*p1ak[j]*p1al[k];
				p3ajkl[i][j][k]=p1aj[i]*p1ak[j]*p1al[k];
			}

	/*initial variables*/
	pairmu[0]=0;	
	pairmu[1]=-elemu[0];
	pairmu[2]=-elemu[1];
	
	paireE1[0][1]=pairen[0];
	paireE1[1][0]=pairen[0];
	paireE1[0][2]=pairen[4];
	paireE1[2][0]=pairen[4];
	paireE1[1][2]=pairen[8];
	paireE1[2][1]=pairen[8];
	
	paireE2[0][1]=pairen[1];
	paireE2[1][0]=pairen[1];
	paireE2[0][2]=pairen[5];
	paireE2[2][0]=pairen[5];
	paireE2[1][2]=pairen[9];
	paireE2[2][1]=pairen[9];
	
	manybodyE[0][1][0][1]=pairen[2];
	manybodyE[1][0][0][1]=pairen[2];
	manybodyE[0][1][1][0]=pairen[2];
	manybodyE[1][0][1][0]=pairen[2];
	
	manybodyE[0][2][0][2]=pairen[6];
	manybodyE[2][0][0][2]=pairen[6];
	manybodyE[0][2][2][0]=pairen[6];
	manybodyE[2][0][2][0]=pairen[6];
	
	manybodyE[1][2][1][2]=pairen[10];
	manybodyE[2][1][1][2]=pairen[10];
	manybodyE[1][2][2][1]=pairen[10];
	manybodyE[2][1][2][1]=pairen[10];
		
	manybodyE[0][1][1][1]=pairen[3];
	manybodyE[1][0][1][1]=pairen[3];
	manybodyE[1][1][0][1]=pairen[3];
	manybodyE[1][1][1][0]=pairen[3];		
	manybodyE[0][2][2][2]=pairen[7];
	manybodyE[2][0][2][2]=pairen[7];
	manybodyE[2][2][0][2]=pairen[7];
	manybodyE[2][2][2][0]=pairen[7];
	manybodyE[1][2][2][2]=pairen[11];
	manybodyE[2][1][2][2]=pairen[11];
	manybodyE[2][2][1][2]=pairen[11];
	manybodyE[2][2][2][1]=pairen[11];

//	fp=fopen("genergyAB0.dat","w");
	
	/*natural iteration */
    while( delta > tau && F < ITERATIONMAX )
	{
		/*vaibles*/
		double lambda;
		double f1x,f2x_a,f2x_b,f3x;
		double p1x,p2x_a,p2x_b,p3x;

			
		/*sum over exp(1/2 lambda beta)*/
		lambda=0.0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						
						
						mu=pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l];
						epsilon=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
						expe=-beta*epsilon+beta*mu/24.0;
					
						if(expe>400)
							energyx=5e+173;
						else 
						{
							if(expe<-pow(10,5))
								energyx=0;
							else 
								energyx=exp(expe);							
						}
						
						f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
						f2x_a=p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l];
						f2x_b=p2aij[i][j]*p2akl[k][l];    
						f3x=p3aijk[i][j][k]*p3aijl[i][j][l]*p3aikl[i][k][l]*p3ajkl[j][k][l];
						
						p1x=pow(f1x,1.0/24);
						p3x=pow(f3x,0.5);
						
						if(f2x_a<1e-70)
							p2x_a=4.641589e+11;
						else 
							p2x_a=pow(f2x_a,-1.0/6);
							
						if(f2x_b<1e-70)
							p2x_b=3.162278e+17;
						else 
							p2x_b=pow(f2x_b,-1.0/4);
						
						lambda+=energyx*p3x*p2x_a*p2x_b*p1x;
							
					}
				}
			}
		}
		
		/*iteration process*/
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						mu=pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l];
						epsilon=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
						expe=-beta*epsilon+beta*mu/24.0;
					
						if(expe>400)
							energyx=5e+173;
						else 
						{
							if(expe<-pow(10,5))
								energyx=0;
							else 
								energyx=exp(expe);							
						}
						
						f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
						f2x_a=p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l];
						f2x_b=p2aij[i][j]*p2akl[k][l];    
						f3x=p3aijk[i][j][k]*p3aijl[i][j][l]*p3aikl[i][k][l]*p3ajkl[j][k][l];
						
						p1x=pow(f1x,1.0/24);
						p3x=pow(f3x,0.5);
						
						if(f2x_a<1e-70)
							p2x_a=4.641589e+11;
						else 
							p2x_a=pow(f2x_a,-1.0/6);
							
						if(f2x_b<1e-70)
							p2x_b=3.162278e+17;
						else 
							p2x_b=pow(f2x_b,-1.0/4);
												
						p4b[i][j][k][l]=energyx*p3x*p2x_a*p2x_b*p1x/lambda;
							
					}
				}
			}
		}
		
		/*generate 3 point site probabilities*/
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
				{
					/*initial*/
					p3bijk[i][j][k]=0;
					p3bijl[i][j][k]=0;
					p3bikl[i][j][k]=0;
					p3bjkl[i][j][k]=0;
					for(l=0;l<3;l++)
					{
						p3bijk[i][j][k]+=p4b[i][j][k][l];
						p3bijl[i][j][k]+=p4b[i][j][l][k];
						p3bikl[i][j][k]+=p4b[i][l][j][k];
						p3bjkl[i][j][k]+=p4b[l][i][j][k];
					}
				}
		
		/*initial p3b,which is average 3 probabilities among 6 sites */
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					p3b[i][j][k]=(p3bijk[i][j][k]+p3bijl[i][j][k]+p3bikl[i][j][k]+p3bjkl[i][j][k])/4.0;
		
		/*generate 2 point site probablilities*/
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
			{	
				/*initial*/
				p2bij[i][j]=0;
				p2bik[i][j]=0;
				p2bil[i][j]=0;
				p2bjk[i][j]=0;
				p2bjl[i][j]=0;
				p2bkl[i][j]=0;
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
					{
						p2bij[i][j]+=p4b[i][j][k][l];
						p2bik[i][j]+=p4b[i][k][j][l];
						p2bil[i][j]+=p4b[i][k][l][j];
						p2bjk[i][j]+=p4b[k][i][j][l];
						p2bjl[i][j]+=p4b[k][i][l][j];
						p2bkl[i][j]+=p4b[k][l][i][j];
					}
			}
		
		/*initial p2b, which is average paire probabilities among 6 sites */
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
			{
				p2b[i][j]=(p2bij[i][j]+p2bik[i][j]+p2bil[i][j]+p2bjk[i][j]+p2bjl[i][j]+p2bkl[i][j])/6.0;
			}
		
		/*generate 1 point site probabliliries*/
		for(i=0;i<3;i++)
		{
			p1bi[i]=0;
			p1bj[i]=0;
			p1bk[i]=0;
			p1bl[i]=0;
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
					{
						p1bi[i]+=p4b[i][j][k][l];
						p1bj[i]+=p4b[j][i][k][l];
						p1bk[i]+=p4b[j][k][i][l];
						p1bl[i]+=p4b[j][k][l][i];
					}
		}
		/*initial p1b*/
		for(i=0;i<3;i++)
			p1b[i]=(p1bi[i]+p1bj[i]+p1bk[i]+p1bl[i])/4.0;
		
		/*the number of iteration */
		F=F+1;
		
		
		/*give the difference of two iterations*/
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				deltaP[3*i+j]=p2a[i][j]-p2b[i][j];
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					deltaP[9+3*(3*i+j)+k]=p3a[i][j][k]-p3b[i][j][k];
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
						deltaP[27+3*(3*(3*i+j)+k)+l]=p4a[i][j][k][l]-p4b[i][j][k][l];
		
		/*caclulate the absnorm deltap*/
		delta=0;
		for(i=0;i<117;i++)
			delta+=fabs(deltaP[i]);
		
	//	printf("%e\n",delta);
		/*reassignment the probabilities, for next loop*/
		for(i=0;i<3;i++)
		{	
			//point probablilites
		    p1ai[i]=p1bi[i];
		    p1aj[i]=p1bj[i];
		    p1ak[i]=p1bk[i];
		    p1al[i]=p1bl[i];
		    //average point probabilities 
		    p1a[i]=p1b[i];
		}
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
		    {	
		    	//pair probabilities  
				p2aij[i][j]=p2bij[i][j];
				p2aik[i][j]=p2bik[i][j];
				p2ail[i][j]=p2bil[i][j];
				p2ajk[i][j]=p2bjk[i][j];
				p2ajl[i][j]=p2bjl[i][j];
				p2akl[i][j]=p2bkl[i][j];
				//average pair probabilities 
				p2a[i][j]=p2b[i][j];
			}
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
				{
					p3aijk[i][j][k]=p3bijk[i][j][k];
					p3aijl[i][j][k]=p3bijl[i][j][k];
					p3aikl[i][j][k]=p3bikl[i][j][k];
					p3ajkl[i][j][k]=p3bjkl[i][j][k];
					//arverage
					p3a[i][j][k]=p3b[i][j][k];
				}
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
						p4a[i][j][k][l]=p4b[i][j][k][l];	
		
	}

	/*grandpotential*/
	double Grandenergy;
	double ge=0;
	double U=0;
	double mutal=0;
	/*internal energy*/
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
				{
					epsilon=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
					U+=6*epsilon*p4b[i][j][k][l];			
				}
				
	/*chemical term*/
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
				{
					mutal+=0.25*p4b[i][j][k][l]*(pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l]);
				}
	/*entrop*/			
	const double Logerr=1e-150;
	ge=0;
	//1 point 
	for(i=0;i<3;i++)
	{
		if(p1bi[i]>Logerr)
		{
			if(p1bj[i]>Logerr)
			{
				if(p1bk[i]>Logerr)
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bj[i]*log(p1bj[i])+p1bk[i]*log(p1bk[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bj[i]*log(p1bj[i])+p1bk[i]*log(p1bk[i]));
				}
				else//k
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bj[i]*log(p1bj[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bj[i]*log(p1bj[i]));
				}
			}
			else//j
			{
				if(p1bk[i]>Logerr)
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bk[i]*log(p1bk[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bk[i]*log(p1bk[i]));
				}
				else//k
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bi[i]*log(p1bi[i]));
				}
			}
		}
		else //i
		{
			if(p1bj[i]>Logerr)
			{
				if(p1bk[i]>Logerr)
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bj[i]*log(p1bj[i])+p1bk[i]*log(p1bk[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bj[i]*log(p1bj[i])+p1bk[i]*log(p1bk[i]));
				}
				else//k
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bj[i]*log(p1bj[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bj[i]*log(p1bj[i]));
				}
			}
			else//j
			{
				if(p1bk[i]>Logerr)
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bk[i]*log(p1bk[i])+p1bl[i]*log(p1bl[i]));
					else				
						ge+=(-0.25)*(p1bk[i]*log(p1bk[i]));
				}
				else//k
				{
					if(p1bl[i]>Logerr)
						ge+=(-0.25)*(p1bl[i]*log(p1bl[i]));
					else				
						ge+=0;
				}
			}
		}
	}
	//2 next nearest neighbor
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			if(p2bij[i][j]>Logerr)
			{
				if(p2bkl[i][j]>Logerr)
					ge+=1.5*(p2bij[i][j]*log(p2bij[i][j])+p2bkl[i][j]*log(p2bkl[i][j]));
				else
					ge+=1.5*(p2bij[i][j]*log(p2bij[i][j]));				
			}
			else
			{
				if(p2bkl[i][j]>Logerr)
					ge+=1.5*(p2bkl[i][j]*log(p2bkl[i][j]));
				else
					ge+=0;	
			}
		}
	//2 Nearest neighbor
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			if(p2bik[i][j]>Logerr)
			{
				if(p2bil[i][j]>Logerr)
				{
					if(p2bjk[i][j]>Logerr)
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bil[i][j]*log(p2bil[i][j])+p2bjk[i][j]*log(p2bjk[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bil[i][j]*log(p2bil[i][j])+p2bjk[i][j]*log(p2bjk[i][j]));
					}
					else//jk
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bil[i][j]*log(p2bil[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bil[i][j]*log(p2bil[i][j]));
					}
				}
				else//il
				{
					if(p2bjk[i][j]>Logerr)
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bjk[i][j]*log(p2bjk[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bjk[i][j]*log(p2bjk[i][j]));
					}
					else//jk
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bik[i][j]*log(p2bik[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bik[i][j]*log(p2bik[i][j]));
					}
				}
			}
			else//ik
			{
				if(p2bil[i][j]>Logerr)
				{
					if(p2bjk[i][j]>Logerr)
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bil[i][j]*log(p2bil[i][j])+p2bjk[i][j]*log(p2bjk[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bil[i][j]*log(p2bil[i][j])+p2bjk[i][j]*log(p2bjk[i][j]));
					}
					else//jk
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bil[i][j]*log(p2bil[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bil[i][j]*log(p2bil[i][j]));
					}
				}
				else//il
				{
					if(p2bjk[i][j]>Logerr)
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bjk[i][j]*log(p2bjk[i][j])+p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=(p2bjk[i][j]*log(p2bjk[i][j]));
					}
					else//jk
					{
						if(p2bjl[i][j]>Logerr)
							ge+=(p2bjl[i][j]*log(p2bjl[i][j]));
						else
							ge+=0;
					}
				}
			}
		}
	//3 point
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
			{
				if(p3bijk[i][j][k]>Logerr)
				{
					if(p3bijl[i][j][k]>Logerr)
					{
						if(p3bikl[i][j][k]>Logerr)
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k]));
						}
						else //ikl
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bijl[i][j][k]*log(p3bijl[i][j][k]));
						}
					}
					else //ijl
					{
						if(p3bikl[i][j][k]>Logerr)
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k]));
						}
						else //ikl
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijk[i][j][k]*log(p3bijk[i][j][k]));
						}
					}
				}
				else //ijk
				{
					if(p3bijl[i][j][k]>Logerr)
					{
						if(p3bikl[i][j][k]>Logerr)
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bikl[i][j][k]*log(p3bikl[i][j][k]));
						}
						else //ikl
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bijl[i][j][k]*log(p3bijl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bijl[i][j][k]*log(p3bijl[i][j][k]));
						}
					}
					else //ijl
					{
						if(p3bikl[i][j][k]>Logerr)
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bikl[i][j][k]*log(p3bikl[i][j][k])+p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=(-3)*(p3bikl[i][j][k]*log(p3bikl[i][j][k]));
						}
						else //ikl
						{
							if(p3bjkl[i][j][k]>Logerr)
								ge+=(-3)*(p3bjkl[i][j][k]*log(p3bjkl[i][j][k]));
							else
								ge+=0;
						}
					}
				}
			}
	//4 point
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
				{	
					if(p4b[i][j][k][l]>Logerr)
						ge+=6*p4b[i][j][k][l]*log(p4b[i][j][k][l]);			
					else 
						ge+=0;
				}
			
	ge=(-1)*KB*ge;
	/*grand energy*/
	Grandenergy=U-tem*ge-mutal;
	
	//returnvalue
	returnvalue[0]=p1bi[0];
	returnvalue[1]=p1bi[1];
	returnvalue[2]=p1bi[2];
	returnvalue[3]=p1bj[0];
	returnvalue[4]=p1bj[1];
	returnvalue[5]=p1bj[2];
	returnvalue[6]=p1bk[0];
	returnvalue[7]=p1bk[1];
	returnvalue[8]=p1bk[2];
	returnvalue[9]=p1bl[0];
	returnvalue[10]=p1bl[1];
	returnvalue[11]=p1bl[2];
	returnvalue[12]=tem;
	returnvalue[13]=elemu[0];
	returnvalue[14]=elemu[1];
	returnvalue[15]=Grandenergy;
	
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
					p4bx[37*i+9*j+3*k+l]=p4b[i][j][k][l];
//	printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f\n",returnvalue[0],returnvalue[1],returnvalue[2],returnvalue[3],returnvalue[4],returnvalue[5],returnvalue[6],returnvalue[7],returnvalue[8],returnvalue[9],returnvalue[10],returnvalue[11],returnvalue[12],returnvalue[13],returnvalue[14],returnvalue[15]);	
}

void bcc3multiproNI( double *paireE,int phase,double tem,double mu1,double mu2,double Tau, double *returnvalue,double *p4bx)
{	
	//give final configuration at fixed tem and mu; multi inipro; 
	//initial probabilites. 
    double p1[NoPROBCONDITION][3]={{0.2,0.4,0.4},{0.4,0.3,0.3},{0.15,0.25,0.6},{0.7,0.2,0.1},{0.3333,0.3333,0.3333}};
    double p2[NoPROBCONDITION][3]={{0.1,0.3,0.6},{0.4,0.2,0.4},{0.1,0.7,0.2},{0.25,0.3,0.45},{0.3,0.4,0.3}};
    double p3[NoPROBCONDITION][3]={{0.3333,0.3333,0.3333},{0.2,0.1,0.7},{0.5,0.25,0.15},{0.3333,0.3333,0.3333},{0.1,0.3,0.6}};
 	//variables 
    double energy[NoPROBCONDITION];
    double NumNI[NoPROBCONDITION];
    double point[NoPROBCONDITION][12];
    double conf[NoPROBCONDITION][16];
    double p4b[NoPROBCONDITION][81];
	double inipro[12];
	double propoint[12];
    int i,j,k;
    
    double mu[2]={mu1,mu2};
    
    for (i=0;i<NoPROBCONDITION;i++)
    {	
    	//initial probabilities 
		if(phase==1)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p1[i][j];
				inipro[6+j]=p1[i][j];
				inipro[9+j]=p1[i][j];
			}
		if(phase==2)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p1[i][j];
				inipro[6+j]=p2[i][j];
				inipro[9+j]=p2[i][j];
			}
		if(phase==3)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p1[i][j];
				inipro[6+j]=p1[i][j];
				inipro[9+j]=p2[i][j];
			}
		
		if(phase==4)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p2[i][j];
				inipro[6+j]=p2[i][j];
				inipro[9+j]=p2[i][j];
			}
		
		if(phase==5)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p2[i][j];
				inipro[6+j]=p1[i][j];
				inipro[9+j]=p2[i][j];
			}
		
		if(phase==6)
			for(j=0;j<3;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[3+j]=p2[i][j];
				inipro[6+j]=p3[i][j];
				inipro[9+j]=p2[NoPROBCONDITION-i][j];
			}

        //NI calculation 
    	bccNI3ele(paireE,inipro,tem,mu,Tau,&conf[i][0],&p4b[i][0]); 
		

    }
  
	int minindex=0;
    for(i=0;i<NoPROBCONDITION;i++)
    {
        if(conf[minindex][15]>=conf[i][15])
            minindex=i;
    }
    
    for(i=0;i<81;i++)
	    p4bx[i]=p4b[minindex][i];
	for(i=0;i<16;i++)
	{
		*(returnvalue+i)=conf[minindex][i];
	}
	for(i=0;i<12;i++)
		propoint[i]=conf[minindex][i];
	*(returnvalue+16)=0.25*(returnvalue[0]+returnvalue[3]+returnvalue[6]+returnvalue[9]);
	*(returnvalue+17)=0.25*(returnvalue[1]+returnvalue[4]+returnvalue[7]+returnvalue[10]);
    *(returnvalue+18)=0.25*(returnvalue[2]+returnvalue[5]+returnvalue[8]+returnvalue[11]); 
	returnvalue[19]=bcc3judgetype(propoint);

}

void bcc3phasetypedata(int s, int region,char data[10] ,char boundary[16], double *paireE, double tem, double *mu1,double *mu2,double Tau,double Granddiff,double Mudiff)
{	
	double index1,index2;
    double conf[5][20];
	double confold[5][20];
	double p4b[5][81];
	int i,j,k;
	double ratio_1,ratio_2;
	struct tm *local;
 	time_t t;
 	 	
	FILE *f1;
	FILE *fp;
	if(s==0)
	{
		f1=fopen(data,"w");
		fp=fopen(boundary,"w");
	}
	else 
	{
		f1=fopen(data,"a");
		fp=fopen(boundary,"a");
	}

	

    	
    int F=0;
    //every tem range
	//all mu range rough scah 
	for(index1=*(mu1+0);index1<=*(mu1+1);index1+=*(mu1+2))
	{	
		if (region==0)
			{
			}
		else
			{
			//insert 
			mu2[0]=2.654*index1+0.57;
			mu2[1]=2.654*index1+1.2;
			mu2[2]=(mu2[1]-mu2[0])/40;
			//insert 
			}
		for(index2=*(mu2+0);index2<=*(mu2+1);index2+=*(mu2+2))
		{	
			if(F%1000==0)
			{
				printf("%d ",F/1000);
				t=time(NULL);
				local=localtime(&t);
    	 		printf(ctime(&t));	
			}
		
			//calculate dis equlilibrum distribution 
			bcc3multiproNI(paireE,1,tem,index1,index2,Tau,&conf[0][0],&p4b[0][0]);
			bcc3multiproNI(paireE,2,tem,index1,index2,Tau,&conf[1][0],&p4b[1][0]);
			bcc3multiproNI(paireE,3,tem,index1,index2,Tau,&conf[2][0],&p4b[2][0]);
			bcc3multiproNI(paireE,4,tem,index1,index2,Tau,&conf[3][0],&p4b[3][0]);
			bcc3multiproNI(paireE,5,tem,index1,index2,Tau,&conf[4][0],&p4b[4][0]);
		 			
		    //sort the three energies
			double energy;
			double point;
			double ni; 
			double rat;
			for(i=0;i<5-1;i++) 
				for(j=i+1;j<5;j++) 
					if(conf[i][15]>conf[j][15]) 
					{ 
						energy=conf[i][15]; 
						conf[i][15]=conf[j][15]; 
						conf[j][15]=energy;
																
						rat=conf[i][16];
						conf[i][16]=conf[j][16];
						conf[j][16]=rat;
						
						rat=conf[i][17];
						conf[i][17]=conf[j][17];
						conf[j][17]=rat;
						
						rat=conf[i][18];
						conf[i][18]=conf[j][18];
						conf[j][18]=rat;
								
						rat=conf[i][19];
						conf[i][19]=conf[j][19];
						conf[j][19]=rat;
						
						for(k=0;k<12;k++)
						{						
							point=conf[i][k];
							conf[i][k]=conf[j][k];
							conf[j][k]=point;				 
						}
						for(k=0;k<81;k++)
						{						
							point=p4b[i][k];
							p4b[i][k]=p4b[j][k];
							p4b[j][k]=point;				 
						}
					} 	
		
			fprintf(f1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][0],conf[0][1],conf[0][2],conf[0][3],conf[0][4],conf[0][5],conf[0][6],conf[0][7],conf[0][8],conf[0][9],conf[0][10],conf[0][11],conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19]);	
	      	//fprintf(f1,"%.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19]);	
	        
			//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][0],p4b[0][1],p4b[0][2],p4b[0][3],p4b[0][4],p4b[0][5],p4b[0][6],p4b[0][7],p4b[0][8]);
	        //fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][9],p4b[0][10],p4b[0][11],p4b[0][12],p4b[0][13],p4b[0][14],p4b[0][15],p4b[0][16],p4b[0][17]);
	       	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][18],p4b[0][19],p4b[0][20],p4b[0][21],p4b[0][22],p4b[0][23],p4b[0][24],p4b[0][25],p4b[0][26]);
	       	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][27],p4b[0][28],p4b[0][29],p4b[0][30],p4b[0][31],p4b[0][32],p4b[0][33],p4b[0][34],p4b[0][35]);
	       	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][36],p4b[0][37],p4b[0][38],p4b[0][39],p4b[0][40],p4b[0][41],p4b[0][42],p4b[0][43],p4b[0][44]);
	    	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][45],p4b[0][46],p4b[0][47],p4b[0][48],p4b[0][79],p4b[0][50],p4b[0][51],p4b[0][52],p4b[0][53]);
	    	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][54],p4b[0][55],p4b[0][56],p4b[0][57],p4b[0][58],p4b[0][59],p4b[0][60],p4b[0][61],p4b[0][62]);
	    	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][63],p4b[0][64],p4b[0][65],p4b[0][66],p4b[0][67],p4b[0][68],p4b[0][69],p4b[0][70],p4b[0][71]);
	    	//fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f ",p4b[0][72],p4b[0][73],p4b[0][74],p4b[0][75],p4b[0][76],p4b[0][77],p4b[0][78],p4b[0][79],p4b[0][80]);	     	
			//fprintf(fp,"%.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19] );

			F=F+1;
			
			
			//restore the result to old one; 
			for(i=0;i<5;i++)
				for(j=0;j<19;j++)
					confold[i][j]=conf[i][j];				
    	}
   			
	}
	
	fclose(fp);
    fclose(f1);

}

int bcc3judgetype(double *p)
{
    int sum=0;
    int i;
    double pi[3],pj[3],pk[3],pl[3];
    double ratio=0;
    const double Errp=2E-4;
    //const double Errp=1E-3;
     for(i=0;i<3;i++)
    {
        pi[i]=p[i];
        pj[i]=p[3+i];
        pk[i]=p[6+i];
        pl[i]=p[9+i];
    }
    for(i=0;i<3;i++)
    {
        if(fabs(pi[i]-pj[i])>Errp)
            sum++;
        if(fabs(pi[i]-pk[i])>Errp)
            sum+=7;
        if(fabs(pi[i]-pl[i])>Errp)
            sum+=7;
        if(fabs(pj[i]-pk[i])>Errp)
            sum+=7;
        if(fabs(pj[i]-pl[i])>Errp)
            sum+=7;
        if(fabs(pk[i]-pl[i])>Errp)
            sum++;
    }
    //printf("%d ",sum);
    switch (sum) {
    case 0:
        return 1;
        break;
    case 30:
        return 3;
        break;
    case 56:
        return 2;
        break;
    case 32:
        return 5;
        break;
    case 58:
        return 7;
        break;
    case 46:
        return 8;
        break;
    case 45:
        return 1;
        break;
    case 72:
        return 7;
        break;
    case 54:
        return 8;
        break;
    case 73:
        return 7;
        break;
    case 61:
        return 8;
        break;
    case 84:
        return 2;
        break;
    case 48:
        return 5;
        break;
    case 86:
        return 7;
        break;
    case 62:
        return 8;
        break;
    case 87:
        return 7;
        break;
    case 69:
        return 8;
        break;
    default:
        return 6;
        break;
    }
}

