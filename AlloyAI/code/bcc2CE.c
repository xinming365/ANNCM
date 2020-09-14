#include <stdio.h>
#include <math.h>
#include <time.h>
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
							8                       0.5
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
void bccNI2ele(double *K,double *G,double *v,double *E,double *pairen, double *inipro, double tem, double elemu, double Tau, double *returnvalue,double *cep);
int phase2elelabel(double *p);
void multiproNI(double *K,double *G,double *v,double *E, double *paireE,double *Elastic_c1,double *Elastic_c2,int phase,double tem,double mu,double Tau, double *returnvalue,double *cep);
void phasetypedata(int s, char data[14] ,char boundary[16],double *K,double *G,double *v,double *E,double *paireE, double *Elastic_c1,double *Elastic_c2,double *tem, double *mu,double Tau,double Granddiff,double Mudiff);
void binaryToManybody(double *w,double *interactions);
int judgetype(double *p);

//the pairen parameter input order is w1AB, w2AB, wAbAB, wABBB.

int main()     
{
	//====================FeAl==========================
	double K[4]={12.309076,-5.268623,1.537342,3.361588};
	double G[4]={9.889214,-6.680118,-5.830203,0.990956};
	double v[4]={-0.006602,0.009334,0.020149,0.002503};
	double E[4]={24.121565,-15.581226,-14.71307,2.771606};
	double Elastic_c1[4]={190.706,79.04035,0.317924,208.3383}; //Fe
	double Elastic_c2[4]={61.42143,32.64593,0.274243,83.19773}; //Al
	
	double paireE[4]={-0.09221762,-0.038328291,0.009263967,0.026809306};
	double tem[3]={100,3700,20};
	double mu[3]={-1.2,1.5,0};
	mu[2]=(mu[1]-mu[0])/186; 
	
	phasetypedata(0,"bccCE10.dat","bccbound10.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,tem,mu,1e-8,1e-5,1e-5);
	//======================================
	return 0;
}


void bccNI2ele(double *K,double *G,double *v,double *E,double *pairen, double *inipro, double tem, double elemu, double Tau, double *returnvalue,double *cep)
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
	double Kfour[16];
	double Gfour[16];
	double vfour[16];
	double Efour[16];
	
	int i,j,k,l;
	double paireE1[2][2]={0};
	double paireE2[2][2]={0};
	double manybodyE[2][2][2][2]={0};
	double pairmu[2]={0};
//	
//	FILE *fp;
//	fp=fopen("bccenergy.dat","w");
	/*declaration of NI variables*/
	int F=ITERATIONMIN+3; 	//the number of iteratin, started from zero. 
	double p3aijk[2][2][2],p3bijk[2][2][2];
	double p3aijl[2][2][2],p3bijl[2][2][2];
	double p3aikl[2][2][2],p3bikl[2][2][2];
	double p3ajkl[2][2][2],p3bjkl[2][2][2];
	double p2aij[2][2],p2bij[2][2];
	double p2aik[2][2],p2bik[2][2];
	double p2ail[2][2],p2bil[2][2];
	double p2ajk[2][2],p2bjk[2][2];
	double p2ajl[2][2],p2bjl[2][2];
	double p2akl[2][2],p2bkl[2][2];
	double p1ai[2],p1bi[2];
	double p1aj[2],p1bj[2];
	double p1ak[2],p1bk[2];
	double p1al[2],p1bl[2];
	
	double p1a[2];
	double p1b[2]={0};
	double p2a[2][2];
	double p2b[2][2]={0};
	double p3a[2][2][2];
	double p3b[2][2][2]={0};
	double p4b[2][2][2][2];
    double p4a[2][2][2][2]={0};
	double deltaP[28];
	double delta=1;
	
	//initial p1a,p2a,p3a.
	for(i=0;i<2;i++)
	{
		p1ai[i]=inipro[i];
		p1aj[i]=inipro[2+i];
		p1ak[i]=inipro[4+i];
		p1al[i]=inipro[6+i];
	}
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
		{
			p2aij[i][j]=p1ai[i]*p1aj[j];
			p2aik[i][j]=p1ai[i]*p1ak[j];
			p2ail[i][j]=p1ai[i]*p1al[j];
			p2ajk[i][j]=p1aj[i]*p1ak[j];
			p2ajl[i][j]=p1aj[i]*p1al[j];
			p2akl[i][j]=p1ak[i]*p1al[j];
	
		}
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
			{
				p3aijk[i][j][k]=p1ai[i]*p1aj[j]*p1ak[k];
				p3aijl[i][j][k]=p1ai[i]*p1aj[j]*p1al[k];
				p3aikl[i][j][k]=p1ai[i]*p1ak[j]*p1al[k];
				p3ajkl[i][j][k]=p1aj[i]*p1ak[j]*p1al[k];
			}
	
	//initial interactions
	binaryToManybody(K,Kfour);
	binaryToManybody(G,Gfour);
	binaryToManybody(v,vfour);
	binaryToManybody(E,Efour);

	/*initial variables*/
	pairmu[0]=elemu;	
	
	paireE1[0][1]=pairen[0];
	paireE1[1][0]=pairen[0];
	
	paireE2[0][1]=pairen[1];
	paireE2[1][0]=pairen[1];
	
	manybodyE[0][1][0][1]=pairen[2];
	manybodyE[1][0][1][0]=pairen[2];
	manybodyE[0][1][1][0]=pairen[2];
	manybodyE[1][0][0][1]=pairen[2];
	
	manybodyE[0][1][1][1]=pairen[3];
	manybodyE[1][0][1][1]=pairen[3];
	manybodyE[1][1][0][1]=pairen[3];
	manybodyE[1][1][1][0]=pairen[3];

	/*natural iteration */
    while( delta > tau && F < ITERATIONMAX && F>ITERATIONMIN)
	{
		/*vaibles*/
		double lambda;
		double f1x,f2x_a,f2x_b,f3x;
		double p1x,p2x_a,p2x_b,p3x;

			
		/*sum over exp(1/2 lambda beta)*/
		lambda=0.0;
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				for(k=0;k<2;k++)
				{
					for(l=0;l<2;l++)
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
		for(i=0;i<2;i++)
		{
			for(j=0;j<2;j++)
			{
				for(k=0;k<2;k++)
				{
					for(l=0;l<2;l++)
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
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
				{
					/*initial*/
					p3bijk[i][j][k]=0;
					p3bijl[i][j][k]=0;
					p3bikl[i][j][k]=0;
					p3bjkl[i][j][k]=0;
					for(l=0;l<2;l++)
					{
						p3bijk[i][j][k]+=p4b[i][j][k][l];
						p3bijl[i][j][k]+=p4b[i][j][l][k];
						p3bikl[i][j][k]+=p4b[i][l][j][k];
						p3bjkl[i][j][k]+=p4b[l][i][j][k];
					}
				}
		
		/*initial p3b,which is average 3 probabilities among 6 sites */
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					p3b[i][j][k]=(p3bijk[i][j][k]+p3bijl[i][j][k]+p3bikl[i][j][k]+p3bjkl[i][j][k])/4.0;
		
		/*generate 2 point site probablilities*/
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
			{	
				/*initial*/
				p2bij[i][j]=0;
				p2bik[i][j]=0;
				p2bil[i][j]=0;
				p2bjk[i][j]=0;
				p2bjl[i][j]=0;
				p2bkl[i][j]=0;
				for(k=0;k<2;k++)
					for(l=0;l<2;l++)
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
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
			{
				p2b[i][j]=(p2bij[i][j]+p2bik[i][j]+p2bil[i][j]+p2bjk[i][j]+p2bjl[i][j]+p2bkl[i][j])/6.0;
			}
		
		/*generate 1 point site probabliliries*/
		for(i=0;i<2;i++)
		{
			p1bi[i]=0;
			p1bj[i]=0;
			p1bk[i]=0;
			p1bl[i]=0;
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					for(l=0;l<2;l++)
					{
						p1bi[i]+=p4b[i][j][k][l];
						p1bj[i]+=p4b[j][i][k][l];
						p1bk[i]+=p4b[j][k][i][l];
						p1bl[i]+=p4b[j][k][l][i];
					}
		}
		/*initial p1b*/
		for(i=0;i<2;i++)
			p1b[i]=(p1bi[i]+p1bj[i]+p1bk[i]+p1bl[i])/4.0;
		
		/*the number of iteration */
		F=F+1;
		
		
		/*give the difference of two iterations*/
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				deltaP[2*i+j]=p2a[i][j]-p2b[i][j];
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					deltaP[4+2*(2*i+j)+k]=p3a[i][j][k]-p3b[i][j][k];
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					for(l=0;l<2;l++)
						deltaP[12+2*(2*(2*i+j)+k)+l]=p4a[i][j][k][l]-p4b[i][j][k][l];
		
		/*caclulate the absnorm deltap*/
		delta=0;
		for(i=0;i<28;i++)
			delta+=fabs(deltaP[i]);
		
	//	printf("%e\n",delta);
		/*reassignment the probabilities, for next loop*/
		for(i=0;i<2;i++)
		{	
			//point probablilites
		    p1ai[i]=p1bi[i];
		    p1aj[i]=p1bj[i];
		    p1ak[i]=p1bk[i];
		    p1al[i]=p1bl[i];
		    //average point probabilities 
		    p1a[i]=p1b[i];
		}
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
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
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
				{
					p3aijk[i][j][k]=p3bijk[i][j][k];
					p3aijl[i][j][k]=p3bijl[i][j][k];
					p3aikl[i][j][k]=p3bikl[i][j][k];
					p3ajkl[i][j][k]=p3bjkl[i][j][k];
					//arverage
					p3a[i][j][k]=p3b[i][j][k];
				}
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					for(l=0;l<2;l++)
						p4a[i][j][k][l]=p4b[i][j][k][l];	
		
	}

	/*grandpotential*/
	double Grandenergy;
	double ge=0;
	double U=0;
	double mutal=0;
	/*internal energy*/
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
				for(l=0;l<2;l++)
				{
					epsilon=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
					U+=6*epsilon*p4b[i][j][k][l];			
				}
				
	/*chemical term*/
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
				for(l=0;l<2;l++)
				{
					mutal+=0.25*p4b[i][j][k][l]*(pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l]);
				}
	/*entrop*/			
	const double Logerr=1e-150;
	ge=0;
	//1 point 
	for(i=0;i<2;i++)
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
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
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
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
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
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
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
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
				for(l=0;l<2;l++)
				{	
					if(p4b[i][j][k][l]>Logerr)
						ge+=6*p4b[i][j][k][l]*log(p4b[i][j][k][l]);			
					else 
						ge+=0;
				}

	ge=(-1)*KB*ge;
	/*grand energy*/
	Grandenergy=U-tem*ge-mutal;
	
	for(i=0;i<4;i++)
		cep[i]=0;
					
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
				for(l=0;l<2;l++)
					{
					cep[0]+=6*Kfour[8*i+4*j+2*k+l]*p4b[i][j][k][l];
					cep[1]+=6*Gfour[8*i+4*j+2*k+l]*p4b[i][j][k][l];
					cep[2]+=6*vfour[8*i+4*j+2*k+l]*p4b[i][j][k][l];
					cep[3]+=6*Efour[8*i+4*j+2*k+l]*p4b[i][j][k][l];
					}
	
	//returnvalue
	returnvalue[0]=p1bi[0];
	returnvalue[1]=p1bi[1];
	returnvalue[2]=p1bj[0];
	returnvalue[3]=p1bj[1];
	returnvalue[4]=p1bk[0];
	returnvalue[5]=p1bk[1];
	returnvalue[6]=p1bl[0];
	returnvalue[7]=p1bl[1];
	returnvalue[8]=tem;
	returnvalue[9]=elemu;
	returnvalue[10]=F;
	returnvalue[11]=Grandenergy;
}

int bccphaselabel(double *p)
{	
    int phase_type=0;
    double p1bi[2],p1bj[2],p1bk[2],p1bl[2];
    int i;
    //const double Errp=0.006;
    const double Errp=0.0005;
    //transform the variables
    for(i=0;i<2;i++)
    {
        p1bi[i]=*(p+0);
        p1bj[i]=*(p+2);
        p1bk[i]=*(p+4);
        p1bl[i]=*(p+6);
    }
    //label the phases 
    if((fabs(p1bi[0]-p1bj[0])<Errp)&&(fabs(p1bi[1]-p1bj[1])<Errp))
    {
    	if((fabs(p1bj[0]-p1bk[0])<Errp)&&( fabs(p1bj[1]-p1bk[1])<Errp))
    	{
			if((fabs(p1bk[0]-p1bl[0])<Errp)&&( fabs(p1bk[1]-p1bl[1])<Errp))
			{
				
				if((fabs(p1bi[0]-0.5)<Errp)&&(fabs(p1bi[1]-0.5)<Errp))
					phase_type=7;
				else
					phase_type=1;
			
			}
			else
			{
				phase_type=3;
			}
		}
		else//i==j,j!=k
		{
			if((fabs(p1bk[0]-p1bl[0])<Errp)&&(fabs(p1bk[1]-p1bl[1])<Errp))
			{
				phase_type=2;
			}
			else
			{
				if((fabs(p1bl[0]-p1bi[0])<Errp)&& (fabs(p1bl[1]-p1bi[1])<Errp))
				{
					phase_type=3;
				}
				else
				{
					phase_type=6;
				}
			}
		}
	}
    else //i!=j
    {
    	if((fabs(p1bj[0]-p1bk[0])<Errp)&&( fabs(p1bj[1]-p1bk[1])<Errp))
    	{
			if((fabs(p1bk[0]-p1bl[0])<Errp)&&( fabs(p1bk[1]-p1bl[1])<Errp))
			{
				phase_type=4;
			}
			else
			{
				if((fabs(p1bl[0]-p1bi[0])<Errp)&& (fabs(p1bl[1]-p1bi[1])<Errp))
				{
					phase_type=5;
				}
				else
				{
					phase_type=9;
				}
			}
		}
		else//i!=j,j!=k
		{
			if((fabs(p1bk[0]-p1bl[0])<Errp)&&(fabs(p1bk[1]-p1bl[1])<Errp))
			{
				if((fabs(p1bl[0]-p1bi[0])<Errp)&& (fabs(p1bl[1]-p1bi[1])<Errp))
				{
					phase_type=4;
				}
				else
				{
					phase_type=9;
				}
			}
			else
			{
				if((fabs(p1bl[0]-p1bi[0])<Errp)&& (fabs(p1bl[1]-p1bi[1])<Errp))
				{
					phase_type=9;
				}
				else
				{
					phase_type=6;
				}
			}
		}
	}    
	
    return phase_type;
}

void multiproNI(double *K,double *G,double *v,double *E, double *paireE,double *Elastic_c1,double *Elastic_c2,int phase,double tem,double mu,double Tau, double *returnvalue,double *cep)
{	
	//give final configuration at fixed tem and mu; multi inipro; 
	//initial probabilites. 
    double p1[NoPROBCONDITION][2]={{0.2,0.8},{0.8,0.2},{0.6,0.4},{0.4,0.6},{0.5,0.5}};
    double p2[NoPROBCONDITION][2]={{0.8,0.2},{0.4,0.6},{0.1,0.9},{0.7,0.3},{0.6,0.4}};
    double p3[NoPROBCONDITION][2]={{0.6,0.4},{0.5,0.5},{0.3,0.7},{0.5,0.5},{0.8,0.2}};
 	//variables 
    double energy[NoPROBCONDITION];
    double NumNI[NoPROBCONDITION];
    double point[NoPROBCONDITION][8];
    double conf[NoPROBCONDITION][12];
    double cepx[NoPROBCONDITION][4];
	double inipro[8];
	double propoint[8];
    int i,j,k;

    for (i=0;i<NoPROBCONDITION;i++)
    {	
    	//initial probabilities 
		if(phase==1)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p1[i][j];
				inipro[4+j]=p1[i][j];
				inipro[6+j]=p1[i][j];
			}
		if(phase==2)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p1[i][j];
				inipro[4+j]=p2[i][j];
				inipro[6+j]=p2[i][j];
			}
		if(phase==3)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p1[i][j];
				inipro[4+j]=p1[i][j];
				inipro[6+j]=p2[i][j];
			}
		
		if(phase==4)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p2[i][j];
				inipro[4+j]=p2[i][j];
				inipro[6+j]=p2[i][j];
			}
		
		if(phase==5)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p2[i][j];
				inipro[4+j]=p1[i][j];
				inipro[6+j]=p2[i][j];
			}
		
		if(phase==6)
			for(j=0;j<2;j++)
			{
				inipro[0+j]=p1[i][j];
				inipro[2+j]=p2[i][j];
				inipro[4+j]=p3[i][j];
				inipro[6+j]=p2[NoPROBCONDITION-i][j];
			}

        //NI calculation 
    	bccNI2ele(K,G,v,E,paireE,inipro,tem,mu,Tau,&conf[i][0],&cepx[i][0]); 
		
    }
	int minindex=0;
    for(i=0;i<NoPROBCONDITION;i++)
    {
        if(conf[minindex][11]>=conf[i][11])
            minindex=i;
    }
    
    
	for(i=0;i<12;i++)
	{
		*(returnvalue+i)=conf[minindex][i];
	}
	for(i=0;i<8;i++)
		propoint[i]=conf[minindex][i];
	*(returnvalue+12)=0.25*(returnvalue[0]+returnvalue[2]+returnvalue[4]+returnvalue[6]);
//	*(returnvalue+13)=bccphaselabel(propoint);       
	returnvalue[13]=judgetype(propoint);

	/*properties*/
    for(i=0;i<4;i++)
	    cep[i]=cepx[minindex][i]+returnvalue[12]*Elastic_c1[i]+(1-returnvalue[12])*Elastic_c2[i];     

}

void binaryToManybody(double *w,double *interactions)
{

	double paireE1[2][2]={0};
	double paireE2[2][2]={0};
	double manybodyE[2][2][2][2]={0};
	int i,j,k,l;
	
	
	paireE1[0][1]=w[0];
	paireE1[1][0]=w[0];
	
	paireE2[0][1]=w[1];
	paireE2[1][0]=w[1];
	
	manybodyE[0][1][0][1]=w[2];
	manybodyE[1][0][1][0]=w[2];
	manybodyE[0][1][1][0]=w[2];
	manybodyE[1][0][0][1]=w[2];
	
	manybodyE[0][1][1][1]=w[3];
	manybodyE[1][0][1][1]=w[3];
	manybodyE[1][1][0][1]=w[3];
	manybodyE[1][1][1][0]=w[3];
	
	for(i=0;i<2;i++)
		for(j=0;j<2;j++)
			for(k=0;k<2;k++)
				for(l=0;l<2;l++)
					interactions[8*i+4*j+2*k+l]=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
					
}


void phasetypedata(int s, char data[14] ,char boundary[16], double *K,double *G,double *v,double *E, double *paireE, double *Elastic_c1,double *Elastic_c2,double *tem, double *mu,double Tau,double Granddiff,double Mudiff)
{	
	double index1,index2;
    double conf[5][14];
	double confold[5][14];
	double cep[5][4]={0};
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

	for (index1=*(tem+0);index1<*(tem+1);index1+=*(tem+2))
    {	
    	if(((int)index1)%10==0)
    	{
			t=time(NULL);
			local=localtime(&t);
    		printf("%d ",(int)index1);
    		printf(ctime(&t));
		}
		//every tem range
		//all mu range rough scah 
		for(index2=*(mu+0);index2<*(mu+1);index2+=*(mu+2))
		{	
			
			//calculate dis equlilibrum distribution 
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,1,index1,index2,Tau,&conf[0][0],&cep[0][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,2,index1,index2,Tau,&conf[1][0],&cep[1][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,3,index1,index2,Tau,&conf[2][0],&cep[2][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,4,index1,index2,Tau,&conf[3][0],&cep[3][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,5,index1,index2,Tau,&conf[4][0],&cep[4][0]);		
		 			
		    //sort the three energies
			double energy;
			double point;
			double ni; 
			double rat;
			for(i=0;i<4;i++) 
				for(j=i+1;j<5;j++) 
					if(conf[i][11]>conf[j][11]) 
					{ 
						energy=conf[i][11]; 
						conf[i][11]=conf[j][11]; 
						conf[j][11]=energy;
											
						ni=conf[i][10];
						conf[i][10]=conf[j][10];
						conf[j][10]=ni;
						
						rat=conf[i][12];
						conf[i][12]=conf[j][12];
						conf[j][12]=rat;
						
						rat=conf[i][13];
						conf[i][13]=conf[j][13];
						conf[j][13]=rat;
						
						for(k=0;k<8;k++)
						{						
							point=conf[i][k];
							conf[i][k]=conf[j][k];
							conf[j][k]=point;				 
						}
						for(k=0;k<4;k++)
						{
							point=cep[i][k];
							cep[i][k]=cep[j][k];
							cep[j][k]=point;
						}
					} 	
			//if(fabs(conf[0][13]-8)>0.03)
				fprintf(f1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.6f %.0f %.6f %.6f %.6f %.6f\n",conf[0][0],conf[0][1],conf[0][2],conf[0][3],conf[0][4],conf[0][5],conf[0][6],conf[0][7],conf[0][8],conf[0][9],conf[0][12],conf[0][13],cep[0][0],cep[0][1],cep[0][2],cep[0][3]);	
			
		//	else
		//		fprintf(f1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.0f %.6f %.6f %.0f\n",conf[1][0],conf[1][1],conf[1][2],conf[1][3],conf[1][4],conf[1][5],conf[1][6],conf[1][7],conf[1][8],conf[1][9],conf[1][10],conf[1][11],conf[1][12],conf[1][13]);	
			/*
			if(fabs(conf[0][13]-7)<0.03)
			{
				conf[0][13]=2;
			}
			else
			{ 
				if((fabs(conf[0][13]-4)<0.03)||(fabs(conf[0][13]-9)<0.03)||(fabs(conf[0][13]-6)<0.03))
				{	
					conf[0][13]=3;
				}
				else 
				{ }
			}
			*/
			//two neibor mu check the lowest energy label,for all ;
			if(index2>*(mu+0))
			{	
				if(fabs(confold[0][13]-conf[0][13])<0.03)
				{
					//two phase is equal, not phase boundary,nothing to do 
				}
				else if(fabs(confold[0][13]-conf[0][13])>0.99)
				{	
					//two phase not equal.  phase boundary
					fprintf(fp,"%.0f %.3f %.3f %.3f %.0f %.0f\n",conf[0][8],conf[0][9],confold[0][12],conf[0][12],confold[0][13],conf[0][13]);	
				}
			}
			
			//restore the result to old one; 
			for(i=0;i<5;i++)
				for(j=0;j<14;j++)
					confold[i][j]=conf[i][j];	
			
    	}

	}
	fclose(fp);
    fclose(f1);
}

int judgetype(double *p)
{
    int sum=0;
    double pi[2],pj[2],pk[2],pl[2];
    int i;
    double ratio=0;
    //const double Errp=0.006;
    const double Errp=3E-3;
    //transform the variables
    for(i=0;i<4;i++)
    	ratio+=*(p+2*i);
    
	for(i=0;i<2;i++)
    {
        pi[i]=p[i];
        pj[i]=p[2+i];
        pk[i]=p[4+i];
        pl[i]=p[6+i];
    }
    
    if(fabs(ratio-0.5)>0.005)
    {
	for(i=0;i<2;i++)
    {
        if(fabs(pi[i]-pj[i])>Errp)
            sum+=2;
        if(fabs(pi[i]-pk[i])>Errp)
            sum+=3;
        if(fabs(pi[i]-pl[i])>Errp)
            sum+=3;
        if(fabs(pj[i]-pk[i])>Errp)
            sum+=3;
        if(fabs(pj[i]-pl[i])>Errp)
            sum+=3;
        if(fabs(pk[i]-pl[i])>Errp)
            sum+=2;
    }
    switch (sum){
    case 0:
        return 1;
        break;
    case 16:
        return 3;
        break;
    case 24:
        return 2;
        break;
    case 20:
        return 5;
        break;
    case 28:
        return 7;
        break;
    case 32:
    	return 6;
    	break; 
    default:
        return 0;
        break;
    }
	}
	else
	{
		return 8;
	}
    
}



