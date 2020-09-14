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
void bccNI3ele(double *K,double *G,double *v,double *E,double *pairen, double *inipro, double tem, double *elemu, double Tau, double *returnvalue,double *cep,double *pair_probability);
void multiproNI(double *K,double *G,double *v,double *E, double *paireE,double *Elastic_c1,double *Elastic_c2,double *Elastic_c3,int phase,double tem,double mu1,double mu2,double Tau, double *returnvalue,double *cep,double *pair_probability);
void phasetypedata(int s,int region, char data[10] ,char boundary[16],double *K,double *G,double *v,double *E,double *paireE, double *Elastic_c1,double *Elastic_c2,double *Elastic_c3,double tem, double *mu1,double *mu2,double Tau,double Granddiff,double Mudiff);
void ThreebinaryToManybody(double propteries[12],double interactions[81]);
void range(double min,double max, int number, double *rangef);
int judgetype(double *p);

//the pairen parameter input order is w1AB, w2AB, wAbAB, wABBB, w1AC,w2AC,.. W1BC, w2BC,, ...
//   ABC is FeCrAl
int main()     
{
	double paireE[12]={0.048040068,-0.023058161,0.010086594,0.009186436,-0.09221762,-0.038328291,0.009263967,0.026809306,-0.018562867,0.00035774,-0.013367108,0.018484173};
	double K[12]={-7.424138,-1.533161,4.375201,0.797415,12.309076,-5.268623,1.537342,3.361588,-8.181009,10.017206,-0.707314,-4.77616};
	double G[12]={-13.079956,19.437046,-14.795175,0.306971,9.889214,-6.680118,-5.830203,0.990956,4.685180,-25.303533,12.819605,16.007890};
	double v[12]={0.021277,-0.041439,0.031523,0.000607,-0.006602,0.009334,0.020149,0.002503,-0.018875,0.066794,-0.026953,-0.056285};
	double E[12]={-32.895939,46.783783,-36.111901,0.844869,24.121565,-15.581226,-14.71307,2.771606,8.433526,-57.298269,30.059156,29.39281};
	//four elasitc properties of 3 constituents. 
	double Elastic_c1[4]={190.706,79.04035,0.317924,208.3383}; //Fe
	double Elastic_c2[4]={261.4871,104.3993,0.323821,276.4119}; //Cr
	double Elastic_c3[4]={61.42143,32.64593,0.274243,83.19773}; //Al
	double T;
	double mu1[3],mu2[3];
	int N=89;
	double Kfour[81];
	int i;
	
	//===========T=600K=======================
	range(-0.75,0.3,N,mu1);
	range(-1.05,1.26,N,mu2);
	phasetypedata(0,0,"bcczfe3a.dat","bcczf3a.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,700,mu1,mu2,1e-8,1e-5,1e-5);
	printf("finished 600K\n");
	
	//===============T=800,1000===========================
	range(-0.92,0.53,N,mu1);
	range(-1.05,1.5,N,mu2);
    phasetypedata(0,0,"bcczfe3b.dat","bcczf3b.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,900,mu1,mu2,1e-8,1e-5,1e-5);
    phasetypedata(0,0,"bcczfe3c.dat","bcczf3c.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,1100,mu1,mu2,1e-8,1e-5,1e-5);
    printf("finished 1000K\n");
    
    
    //================T=1400====================================
    range(-1,0.5,N,mu1);
    range(-1.1,1.65,N,mu2);
            
    phasetypedata(0,0,"bcczfe3d.dat","bcczf3d.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,1200,mu1,mu2,1e-8,1e-5,1e-5);
    printf("finished 1400K\n ");
	
	//===============T=2000=====================================
//	range(-1.3,0.8,N,mu1);
//	range(-1.2,1.8,N,mu2);
//    phasetypedata(0,0,"bcczfe3z.dat","bcczf3z.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,2000,mu1,mu2,1e-8,1e-5,1e-5);
//    printf("finished 2000K\n");	
	
	/*
	//===========T=3000===========================================
	range(-1.4,1.2,N,mu1);
	range(-1.6,2.3,N,mu2);
	phasetypedata(0,0,"bcczx.dat","bcczf3x.dat",K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,3000,mu1,mu2,1e-8,1e-5,1e-5);
    printf("finished 3000K\n");	
	*/
	//sound
	Beep(784, 700);
	Beep(880, 800);
	Beep(932, 900);
	Beep(1046, 1000);
	Beep(1175, 1100);
	
	return 0;
}

void range(double min,double max, int number, double *rangef)
{
	rangef[0]=min;
	rangef[1]=max;
	rangef[2]=(max-min)/number;
}

void bccNI3ele(double *K,double *G,double *v,double *E, double *pairen, double *inipro, double tem, double *elemu, double Tau, double *returnvalue,double *cep,double *pair_probability)
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
	double Kfour[81];
	double Gfour[81];
	double vfour[81];
	double Efour[81];
	
	
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
	
	
	//initial interactions
	ThreebinaryToManybody(K,Kfour);
	ThreebinaryToManybody(G,Gfour);
	ThreebinaryToManybody(v,vfour);
	ThreebinaryToManybody(E,Efour);

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
	
	for(i=0;i<4;i++)
		cep[i]=0;
	/*excess properties*/
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
				{
					cep[0]+=6*Kfour[27*i+9*j+3*k+l]*p4b[i][j][k][l];
					cep[1]+=6*Gfour[27*i+9*j+3*k+l]*p4b[i][j][k][l];
					cep[2]+=6*vfour[27*i+9*j+3*k+l]*p4b[i][j][k][l];
					cep[3]+=6*Efour[27*i+9*j+3*k+l]*p4b[i][j][k][l];
				}
	
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
	
	pair_probability[0]=p2bij[0][0];
	pair_probability[1]=p2bij[0][1];
	pair_probability[2]=p2bij[0][2];
	pair_probability[3]=p2bij[1][0];
	pair_probability[4]=p2bij[1][1];
	pair_probability[5]=p2bij[1][2];
	pair_probability[6]=p2bij[2][0];
	pair_probability[7]=p2bij[2][1];
	pair_probability[8]=p2bij[2][2];
	
	pair_probability[9]=p2bik[0][0];
	pair_probability[10]=p2bik[0][1];
	pair_probability[11]=p2bik[0][2];
	pair_probability[12]=p2bik[1][0];
	pair_probability[13]=p2bik[1][1];
	pair_probability[14]=p2bik[1][2];
	pair_probability[15]=p2bik[2][0];
	pair_probability[16]=p2bik[2][1];
	pair_probability[17]=p2bik[2][2];
	
	pair_probability[18]=p2bil[0][0];
	pair_probability[19]=p2bil[0][1];
	pair_probability[20]=p2bil[0][2];
	pair_probability[21]=p2bil[1][0];
	pair_probability[22]=p2bil[1][1];
	pair_probability[23]=p2bil[1][2];
	pair_probability[24]=p2bil[2][0];
	pair_probability[25]=p2bil[2][1];
	pair_probability[26]=p2bil[2][2];

	pair_probability[27]=p2bjk[0][0];
	pair_probability[28]=p2bjk[0][1];
	pair_probability[29]=p2bjk[0][2];
	pair_probability[30]=p2bjk[1][0];
	pair_probability[31]=p2bjk[1][1];
	pair_probability[32]=p2bjk[1][2];
	pair_probability[33]=p2bjk[2][0];
	pair_probability[34]=p2bjk[2][1];
	pair_probability[35]=p2bjk[2][2];
	
	pair_probability[36]=p2bjl[0][0];
	pair_probability[37]=p2bjl[0][1];
	pair_probability[38]=p2bjl[0][2];
	pair_probability[39]=p2bjl[1][0];
	pair_probability[40]=p2bjl[1][1];
	pair_probability[41]=p2bjl[1][2];
	pair_probability[42]=p2bjl[2][0];
	pair_probability[43]=p2bjl[2][1];
	pair_probability[44]=p2bjl[2][2];
	
	pair_probability[45]=p2bkl[0][0];
	pair_probability[46]=p2bkl[0][1];
	pair_probability[47]=p2bkl[0][2];
	pair_probability[48]=p2bkl[1][0];
	pair_probability[49]=p2bkl[1][1];
	pair_probability[50]=p2bkl[1][2];
	pair_probability[51]=p2bkl[2][0];
	pair_probability[52]=p2bkl[2][1];
	pair_probability[53]=p2bkl[2][2];


//	printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f\n",returnvalue[0],returnvalue[1],returnvalue[2],returnvalue[3],returnvalue[4],returnvalue[5],returnvalue[6],returnvalue[7],returnvalue[8],returnvalue[9],returnvalue[10],returnvalue[11],returnvalue[12],returnvalue[13],returnvalue[14],returnvalue[15]);	
}

void multiproNI(double *K,double *G,double *v,double *E, double *paireE, double *Elastic_c1,double *Elastic_c2,double *Elastic_c3, int phase,double tem,double mu1,double mu2,double Tau, double *returnvalue,double *cep,double *pair_probability)
{	
	//give final configuration at fixed tem and mu; multi inipro; 
	//initial probabilites. 
    //double p1[NoPROBCONDITION][3]={{0.2,0.4,0.4},{0.4,0.3,0.3},{0.15,0.25,0.6},{0.7,0.2,0.1},{0.3333,0.3333,0.3333}};
    //double p2[NoPROBCONDITION][3]={{0.1,0.3,0.6},{0.4,0.2,0.4},{0.1,0.7,0.2},{0.25,0.3,0.45},{0.3,0.4,0.3}};
    double p1[NoPROBCONDITION][3]={{0.2,0.3,0.5},{0.1,0.2,0.7},{0.3,0.3,0.4},{0.25,0.35,0.4},{0.15,0.15,0.7}};
    double p2[NoPROBCONDITION][3]={{0.6,0.2,0.2},{0.1,0.3,0.6},{0.3,0.4,0.3},{0.25,0.3,0.45},{0.15,0.25,0.6}};

    double p3[NoPROBCONDITION][3]={{0.3333,0.3333,0.3333},{0.2,0.1,0.7},{0.5,0.25,0.15},{0.3333,0.3333,0.3333},{0.1,0.3,0.6}};
 	//variables 
    double energy[NoPROBCONDITION];
    double NumNI[NoPROBCONDITION];
    double point[NoPROBCONDITION][12];
    double conf[NoPROBCONDITION][16];
    double cepx[NoPROBCONDITION][4];
    double pair[NoPROBCONDITION][54];
	double inipro[12];
	double propoint[12];
    int i,j,k;
    
    double mu[2]={mu1,mu2};
    
//	FILE *fp;
//	fp=fopen("test.dat","w");
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
    	bccNI3ele(K,G,v,E,paireE,inipro,tem,mu,Tau,&conf[i][0],&cepx[i][0],&pair[i][0]); 
		
//	//	fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f\n",inipro[0],inipro[3],inipro[6],inipro[9],conf[i][0],conf[i][1],conf[i][2],conf[i][3],conf[i][4],conf[i][5],conf[i][6],conf[i][7],conf[i][8],conf[i][9],conf[i][10],conf[i][11],conf[i][12],conf[i][13],conf[i][14],conf[i][15]);	
//		printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f\n",inipro[0],inipro[3],inipro[6],inipro[9],conf[i][0],conf[i][1],conf[i][2],conf[i][3],conf[i][4],conf[i][5],conf[i][6],conf[i][7],conf[i][8],conf[i][9],conf[i][10],conf[i][11],conf[i][12],conf[i][13],conf[i][14],conf[i][15]);	

    }
  
	int minindex=0;
    for(i=0;i<NoPROBCONDITION;i++)
    {
        if(conf[minindex][15]>=conf[i][15])
            minindex=i;
    }
    
	for(i=0;i<16;i++)
	{
		*(returnvalue+i)=conf[minindex][i];
	}
	for(i=0;i<12;i++)
		propoint[i]=conf[minindex][i];
	*(returnvalue+16)=0.25*(returnvalue[0]+returnvalue[3]+returnvalue[6]+returnvalue[9]);
	*(returnvalue+17)=0.25*(returnvalue[1]+returnvalue[4]+returnvalue[7]+returnvalue[10]);
    *(returnvalue+18)=0.25*(returnvalue[2]+returnvalue[5]+returnvalue[8]+returnvalue[11]);
//	*(returnvalue+13)=bccphaselabel(propoint); 
	returnvalue[19]=judgetype(propoint);
	
	
	/*properties*/
    for(i=0;i<4;i++)
	    cep[i]=cepx[minindex][i]+returnvalue[16]*Elastic_c1[i]+returnvalue[17]*Elastic_c2[i]+returnvalue[18]*Elastic_c3[i];    
	for(i=0;i<54;i++)
		pair_probability[i]=pair[minindex][i]; 
//	printf("%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f\n",returnvalue[0],returnvalue[1],returnvalue[2],returnvalue[3],returnvalue[4],returnvalue[5],returnvalue[6],returnvalue[7],returnvalue[8],returnvalue[9],returnvalue[10],returnvalue[11],returnvalue[12],returnvalue[13],returnvalue[14],returnvalue[15]);	
//	fclose(fp);
}

void ThreebinaryToManybody(double propteries[12],double interactions[81])
{
	double temp[4];
	double w1[4],w2[4],w3[4];
	double paireE1[3][3]={0},paireE2[3][3]={0};
	double manybodyE[3][3][3][3]={0};

	int i,j,k,l;

	double f;	
	//initial the binary interactions. 
	for(i=0;i<4;i++)
	{
		w1[i]=propteries[i+0];
		w2[i]=propteries[i+4];
		w3[i]=propteries[i+8];
	}
	
	//initial the ternary four-body interactions
	paireE1[0][1]=w1[0];
	paireE1[1][0]=w1[0];
	paireE1[0][2]=w2[0];
	paireE1[2][0]=w2[0];
	paireE1[1][2]=w3[0];
	paireE1[2][1]=w3[0];
	
	paireE2[0][1]=w1[1];
	paireE2[1][0]=w1[1];
	paireE2[0][2]=w2[1];
	paireE2[2][0]=w2[1];
	paireE2[1][2]=w3[1];
	paireE2[2][1]=w3[1];
	
	manybodyE[0][1][0][1]=w1[2];
	manybodyE[1][0][0][1]=w1[2];
	manybodyE[0][1][1][0]=w1[2];
	manybodyE[1][0][1][0]=w1[2];
	
	manybodyE[0][2][0][2]=w2[2];
	manybodyE[2][0][0][2]=w2[2];
	manybodyE[0][2][2][0]=w2[2];
	manybodyE[2][0][2][0]=w2[2];
	
	manybodyE[1][2][1][2]=w3[2];
	manybodyE[2][1][1][2]=w3[2];
	manybodyE[1][2][2][1]=w3[2];
	manybodyE[2][1][2][1]=w3[2];
		
	manybodyE[0][1][1][1]=w1[3];
	manybodyE[1][0][1][1]=w1[3];
	manybodyE[1][1][0][1]=w1[3];
	manybodyE[1][1][1][0]=w1[3];		
	manybodyE[0][2][2][2]=w2[3];
	manybodyE[2][0][2][2]=w2[3];
	manybodyE[2][2][0][2]=w2[3];
	manybodyE[2][2][2][0]=w2[3];
	manybodyE[1][2][2][2]=w3[3];
	manybodyE[2][1][2][2]=w3[3];
	manybodyE[2][2][1][2]=w3[3];
	manybodyE[2][2][2][1]=w3[3];
	
	
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			for(k=0;k<3;k++)
				for(l=0;l<3;l++)
					interactions[27*i+9*j+3*k+l]=(paireE1[i][k]+paireE1[i][l]+paireE1[j][k]+paireE1[j][l])/6.0+(paireE2[i][j]+paireE2[k][l])/4.0+manybodyE[i][j][k][l];
	
}

void phasetypedata(int s, int region,char data[10] ,char boundary[16], double *K,double *G,double *v,double *E,double *paireE, double *Elastic_c1,double *Elastic_c2,double *Elastic_c3,double tem, double *mu1,double *mu2,double Tau,double Granddiff,double Mudiff)
{	
	double index1,index2;
    double conf[5][20];
	double confold[5][20];
	double cep[5][4];
	double pair[5][54];
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
	//all mu range rough scan 
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
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,1,tem,index1,index2,Tau,&conf[0][0],&cep[0][0],&pair[0][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,2,tem,index1,index2,Tau,&conf[1][0],&cep[1][0],&pair[1][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,3,tem,index1,index2,Tau,&conf[2][0],&cep[2][0],&pair[2][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,4,tem,index1,index2,Tau,&conf[3][0],&cep[3][0],&pair[3][0]);
			multiproNI(K,G,v,E,paireE,Elastic_c1,Elastic_c2,Elastic_c3,5,tem,index1,index2,Tau,&conf[4][0],&cep[4][0],&pair[4][0]);
		 			
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
						for(k=0;k<4;k++)
						{						
							point=cep[i][k];
							cep[i][k]=cep[j][k];
							cep[j][k]=point;				 
						}
						
						for(k=0;k<54;k++)
						{						
							point=pair[i][k];
							pair[i][k]=pair[j][k];
							pair[j][k]=point;				 
						}
					} 	
		
			//fprintf(f1,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][0],conf[0][1],conf[0][2],conf[0][3],conf[0][4],conf[0][5],conf[0][6],conf[0][7],conf[0][8],conf[0][9],conf[0][10],conf[0][11],conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19]);	
	      //fprintf(f1,"%.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19]);	
	      	     	
	      	for(i=0;i<54;i++)
	      		fprintf(f1,"%.6f ",pair[0][i]);
	      	fprintf(f1, "%.6f %.6f  %.8f %.6f %.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",cep[0][0],cep[0][1],cep[0][2],cep[0][3],conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19] );
	
			fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f  %.6f %.6f %.8f %.6f %.0f %.3f %.3f %.6f %.6f %.6f %.6f %.0f\n",conf[0][0],conf[0][1],conf[0][2],conf[0][3],conf[0][4],conf[0][5],conf[0][6],conf[0][7],conf[0][8],conf[0][9],conf[0][10],conf[0][11],cep[0][0],cep[0][1],cep[0][2],cep[0][3],conf[0][12],conf[0][13],conf[0][14],conf[0][15],conf[0][16],conf[0][17],conf[0][18],conf[0][19] );
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
			F=F+1;
			
			//two neibor mu check the lowest energy label,for all ;
//			if(index1>*(mu1+0))
//			{	
//				if(fabs(confold[0][19]-conf[0][19])<0.03)
//				{
//					//two phase is equal, not phase boundary,nothing to do 
//				}
//				else if(fabs(confold[0][19]-conf[0][19])>0.99)
//				{	
//					//two phase not equal.  phase boundary
//					fprintf(fp,"%.6f %.6f %.6f %.6f %.6f %.6f\n",conf[0][16],conf[0][17],conf[0][18],confold[0][16],confold[0][17],confold[0][18]);		
//				}
//			}
			
			//restore the result to old one; 
			for(i=0;i<5;i++)
				for(j=0;j<19;j++)
					confold[i][j]=conf[i][j];	
			
    	}
   			
	}
	
	fclose(fp);
    fclose(f1);

}


int judgetype(double *p)
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

