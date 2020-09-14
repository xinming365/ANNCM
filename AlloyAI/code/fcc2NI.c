#include <stdio.h>
#include <math.h>
#include <string.h>

#define KB 8.61173305e-5
#define CHEMICALRANGEMAX 5000
#define ITERATIONMAX 2000
#define ITERATIONMIN 30
#define NoPROBCONDITION 5 
#define Z 12
#define TAU 1e-13
#define ERRP 0.1
#define LOGERR 1e-150
#define DIFFDELTA 1e-6

/*there's some problem in memcpy or strcpy function*/
/*struct of 2 ele configurations and 3 ele configurationss*/
struct Configuration2ele
{
	char Atoms[4]; 	        // the species of element
	int NumIt;              // Number of iterations 
	double Tem;             // temperature
	double Mu;              // chemical potential  
	double Grandenergy;     // grandpotential 
	double pointP[8];       // site probablilities 
};

struct Configuration3ele
{
	char Atoms[6];
	int NumIt;
	double Tem;
	double Mu;
	double Grandenergy;
	double pointP[12];
};

struct Diffenergy
{
	double G_energy;
	int phase_1;
	double p_1[8];
	int phase_2;
	double p_2[8];
};

/*declaration of function*/
struct Configuration2ele NI2ele(char *atoms, double *paireE, double inipro[3][4], char *symmetry, double tem, double elemu);   //AT GIVEN temperature and chemical potential, give the guess pairprobablilities to Iterate final probabilites, return struct Configuration wich contain the final configuration and energy information.
struct Configuration3ele NI3ele(char *atoms, double *paireE, double inipro[3][9], char *symmetry, double tem, double elemu);
void print2ele(struct Configuration2ele conf);                                                                                 //print struct Configuration2ele
void print3ele(struct Configuration3ele conf);  
void phasedata2eleAB(char *atoms,double *paireE, double *tem,double *mu);     
void phasedata2eledisorder(char *atoms,double *paireE, double *tem, double *mu);
void phasedata2eleA3B(char *atoms,double *paireE, double *tem,double *mu);   
void test(char *atoms,double *paireE,double tem, double mu);                                                  
int phase2elelabel(double *p);       //label phases (2 elements)
int phase3elelabel(double *p);       //label phases (3 elements)
struct Diffenergy energydifference2ele(char *atoms, double *paireE, char *s_1, char *s_2, double tem,double mu);   //calculate energydifference of symmetry s1 and s2 phases. at given tem and chemical potential. RETURN energydifference and the phaselabels 
void bisectionroot(char *atoms,double *paireE,char *s_1, char *s_2,double *tem,double *mu);    //find root of diffenergy use the method of bisection 

int main()
{
	char *s;
	double pairen[3]={0.0835909,0.01698049,-0.75723212};
	double tem[3]={50,400,10};
	double mu[3]={-0.2,0.2,0.001}; 
	
	//struct Diffenergy diff1;
	//diff1=energydifference2ele("FeCr",pairen,"A3B","AB",200,0.0990);
	//printf("%f",diff1.G_energy);
	//printf(" %d %d \n",diff1.phase_1,diff1.phase_2);
	
	bisectionroot("FeCr",pairen,"A3B","dis",tem,mu);
	bisectionroot("FeCr",pairen,"AB","dis",tem,mu);
	bisectionroot("FeCr",pairen,"A3B","AB",tem,mu);
	
	/*phasedata*/
    //phasedata2eledisorder("FeCr",pairen, tem, mu);
	phasedata2eleAB("FeCr",pairen, tem, mu);
   	//phasedata2eleA3B("FeCr",pairen, tem, mu);
    
	return 0;
}

struct Configuration2ele NI2ele(char *atoms, double *paireE, double inipro[3][4], char *symmetry, double tem, double elemu)
{	
    //AT GIVEN temperature and chemical potential,
	//give the guess pairprobablilities to Iterate final probabilites
	//return struct Configuration wich contain the final configuration and energy information.
	
	//paireE 1d double array
	//inipro[3][4] 3 sets of inipairprobabilities
	
	/*declaration of variables*/
	double beta=1.0/(tem*KB);
	double tau;
	double pairmu[2]={0};
	double pairen[2][2]={{0},{0}};
	double mun;
	int i,j,k,l;
	
	/*declaration of NI variables*/
	int F=0; 							//the number of iteratin, started from zero. 
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
	double mu,energy,expe,energyx;
	
	double p1a[2];
	double p1b[2]={0};
	double p2a[2][2],p2b[2][2];
	double p4b[2][2][2][2];
	double p4a[2][2][2][2]={0};
	double deltaP[20];
	double delta=1;
	double *point,*pair;
	double Grandpotential[ITERATIONMAX]={0};
	
	struct Configuration2ele conf;
	double inipoint[4][2]={{0},{0},{0},{0}};

	/*initial variables*/
	
	
	//terminate conditition 
	if (fabs(elemu-0.0)>=1e-6)
		tau=TAU*0.01/(1000*pow(10,4.2)*fabs(elemu));
	else
		tau=TAU;

	//the parameters 
	pairmu[0]=*(paireE+0)+elemu;
	pairen[0][1]=*(paireE+1);
	pairen[1][0]=*(paireE+1);
	mun=*(paireE+2); //for grand potential calculation;
	
	//For 3 symmetries, use input pairprobablilities inipro to generate point probabilities 
	//initial p2a, which is the average of pair probabilities among 6 sites  
	if(memcmp("A3B",symmetry,3)==0)
	{	
		for(i=0;i<3;i++)
		{	
			//sum along index 1
			inipoint[i][0]=inipro[0][0]+inipro[0][1]; 
			inipoint[i][1]=inipro[0][2]+inipro[0][3];
		}
		
		//sum along index 0			
		inipoint[3][0]=inipro[0][0]+inipro[0][2];
		inipoint[3][1]=inipro[0][1]+inipro[0][3];
		
		//initial p2a
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				p2a[i][j]=(3*inipro[0][2*i+j]+3*inipro[1][2*i+j])/6.0;
	}
	else if (memcmp("AB",symmetry,2)==0)
	{
		
		for(i=0;i<2;i++)
		{
			//sum along index 1
			inipoint[i][0]=inipro[1][0]+inipro[1][1]; 
			inipoint[i][1]=inipro[1][2]+inipro[1][3];
		}
		for(i=2;i<4;i++)
		{
			//sum along index 0
			inipoint[i][0]=inipro[1][0]+inipro[1][2]; 
			inipoint[i][1]=inipro[1][1]+inipro[1][3];
		}
	 
	 	//initial p2a
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				p2a[i][j]=(inipro[0][2*i+j]+4*inipro[1][2*i+j]+inipro[2][2*i+j])/6.0;
	}
	else 
	{
		for(i=0;i<4;i++)
		{
			//sum along index 1
			inipoint[i][0]=inipro[0][0]+inipro[0][1]; 
			inipoint[i][1]=inipro[0][2]+inipro[0][3];
		}
		
		
		//initial p2a
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				p2a[i][j]=inipro[0][2*i+j];
	}
	
	
	/*natural iteration */
//	fp=fopen("2.dat","w");
    while( delta > tau && F < ITERATIONMAX )
	{
		
		/*vaibles*/
		double lambda; 
		double p2x,p1x,f2x,f1x;
		double epoint=0;
		double epair=0;
		double ge=0;
		
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
						energy=0.5*(pairen[i][j]+pairen[i][k]+pairen[i][l]+pairen[j][k]+pairen[j][l]+pairen[k][l]);
						expe=2*beta*energy-1.5*beta*mu;
						
						if(expe>400)
							energyx=5e+173;
						else if(expe<-1e10)
							energyx=1.0+expe;
						else 
							energyx=exp(expe);							

						
						if(F==0)
						{
							if(memcmp("A3B",symmetry,3)==0)
							{
								p2x=pow(inipro[0][2*i+l]*inipro[0][2*j+l]*inipro[0][2*k+l]*inipro[1][2*i+j]*inipro[1][2*i+k]*inipro[1][2*j+k],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);									
							}
							else if(memcmp("AB",symmetry,2)==0)
							{
								p2x=pow(inipro[0][2*i+j]*inipro[1][2*i+k]*inipro[1][2*i+l]*inipro[1][2*j+k]*inipro[1][2*j+l]*inipro[2][2*k+l],0.5);   
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
							else 
							{
								p2x=pow(inipro[0][2*i+j]*inipro[0][2*i+k]*inipro[0][2*i+l]*inipro[0][2*j+k]*inipro[0][2*j+l]*inipro[0][2*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
						}
						else 
						{	
							f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
							f2x=p2aij[i][j]*p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l]*p2akl[k][l];    
							p2x=pow(f2x,0.5);
							if(f1x<1e-70)
								p1x=5.6234132519e+43;
							else 
								p1x=pow(f1x,-0.625);
						}
							
						lambda+=energyx*p2x*p1x;
							
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
						energy=0.5*(pairen[i][j]+pairen[i][k]+pairen[i][l]+pairen[j][k]+pairen[j][l]+pairen[k][l]);
						expe=2*beta*energy-1.5*beta*mu;
						
						if(expe>400)
							energyx=5e+173;
						else if(expe<-1e10)
							energyx=1.0+expe;
						else 
							energyx=exp(expe);							
				
						if(F==0)
						{
							if(memcmp("A3B",symmetry,3)==0)
							{
								p2x=pow(inipro[0][2*i+l]*inipro[0][2*j+l]*inipro[0][2*k+l]*inipro[1][2*i+j]*inipro[1][2*i+k]*inipro[1][2*j+k],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);									
							}
							else if(memcmp("AB",symmetry,2)==0)
							{
								p2x=pow(inipro[0][2*i+j]*inipro[1][2*i+k]*inipro[1][2*i+l]*inipro[1][2*j+k]*inipro[1][2*j+l]*inipro[2][2*k+l],0.5);   
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
							else 
							{
								p2x=pow(inipro[0][2*i+j]*inipro[0][2*i+k]*inipro[0][2*i+l]*inipro[0][2*j+k]*inipro[0][2*j+l]*inipro[0][2*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
						}
						else 
						{	
							f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
							f2x=p2aij[i][j]*p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l]*p2akl[k][l];   
							p2x=pow(f2x,0.5);
							if(f1x<1e-70)
								p1x=5.6234132519e+43;
							else 
								p1x=pow(f1x,-0.625);
						}
							
						p4b[i][j][k][l]=energyx*p2x*p1x/lambda;
							
					}
				}
			}
		}
		
		
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
					for(l=0;l<2;l++)
						deltaP[4+2*(2*(2*i+j)+k)+l]=p4a[i][j][k][l]-p4b[i][j][k][l];
		
		/*caclulate the absnorm deltap*/
		delta=0;
		for(i=0;i<20;i++)
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
					for(l=0;l<2;l++)
						p4a[i][j][k][l]=p4b[i][j][k][l];
	
				
		
		/*grandpotential*/
		for(i=0;i<2;i++)
			{
				epoint+=p1b[i]*pairmu[i];
			}
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				{
					epair+=0.5*p2b[i][j]*pairen[i][j];
				}
		
		/*entrop log of number of macro state*/
		ge=0;
		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
			{
				if(p2bij[i][j]>=LOGERR)
				{
					if(p2bik[i][j]>=LOGERR)
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1));
								}
							}
						}
					}
					else //ik
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1));
								}
							}
						}
					}
				}
				else // ij
				{
					if(p2bik[i][j]>=LOGERR)
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1));
								}
							}
						}
					}
					else //ik
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=0;
								}
							}
						}
					}
				}	
			}
		for(i=0;i<2;i++)
		{
			if(p1bi[i]>=LOGERR)
			{
				if(p1bj[i]>=LOGERR)
				{
					if(p1bk[i]>=LOGERR)
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1));	
					}
					else
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1));
					}
				}
				else
				{
					if(p1bk[i]>LOGERR)
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bk[i]*(log(p1bk[i])-1));
					}
					else
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1));	
					}
				}
			}
			else
			{
				if(p1bj[i]>=LOGERR)
				{
					if(p1bk[i]>=LOGERR)
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1));	
					}
					else
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1));
					}
				}
				else
				{
					if(p1bk[i]>LOGERR)
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bk[i]*(log(p1bk[i])-1));
					}
					else
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*0;	
					}
				}
			}
		}

		for(i=0;i<2;i++)
			for(j=0;j<2;j++)
				for(k=0;k<2;k++)
					for(l=0;l<2;l++)
					{	
						if(p4b[i][j][k][l]>LOGERR)
							ge+=-2*p4b[i][j][k][l]*(log(p4b[i][j][k][l])-1);			
						else 
							ge+=0;
					}
			
		
		Grandpotential[F-1]=Z*epoint-2*Z*epair+Z*mun-KB*tem*ge;
		
		/*begin inserted */
//		int phase;
//		double pp[8];
//		pp[0]=p1bi[0];
//		pp[1]=p1bi[1];
//		pp[2]=p1bj[0];
//		pp[3]=p1bj[1];
//		pp[4]=p1bk[0];
//		pp[5]=p1bk[1];
//		pp[6]=p1bl[0];
//		pp[7]=p1bl[1];
//		phase=phase2elelabel(pp);
//		/*end inserted*/
//		fprintf(fp,"%8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f %.16f %d\n",p1bi[0],p1bi[1],p1bj[0],p1bj[1],p1bk[0],p1bk[1],p1bl[0],p1bl[1],Grandpotential[i],phase);	

	
	}
	
//	fclose(fp);
	
	
	/*initial return value*/	
	memcpy(conf.Atoms,atoms,4);
	conf.NumIt=F;
	conf.Tem=tem;
	conf.Mu=elemu;
	conf.Grandenergy=Grandpotential[F-1];
	conf.pointP[0]=p1bi[0];
	conf.pointP[1]=p1bi[1];
	conf.pointP[2]=p1bj[0];
	conf.pointP[3]=p1bj[1];
	conf.pointP[4]=p1bk[0];
	conf.pointP[5]=p1bk[1];
	conf.pointP[6]=p1bl[0];
	conf.pointP[7]=p1bl[1];	
	
	return conf;
}

void print2ele(struct Configuration2ele conf)
{
	int i;
	printf("Conf Atoms: %s\n",conf.Atoms);
	printf("iteration number is: %d\n",conf.NumIt);
	printf("the temperature is %.2f\n",conf.Tem);
	printf("the chemical potential is %.4f\n",conf.Mu);
	printf("grand potential is %.4f\n",conf.Grandenergy);
	printf("the probablility is:\n");
	for(i=0;i<4;i++)
		printf("%8.3f%8.3f\n",conf.pointP[2*i],conf.pointP[2*i+1]);
	
}

struct Configuration3ele NI3ele(char *atoms, double *paireE, double inipro[3][9], char *symmetry, double tem, double elemu)
{	
    //AT GIVEN temperature and chemical potential,
	//give the guess pairprobablilities to Iterate final probabilites
	//return struct Configuration wich contain the final configuration and energy information.
	
	//paireE 1d double array
	//inipro[3][4] 3 sets of inipairprobabilities
	
	/*declaration of variables*/
	double beta=1/(tem*KB);
	double tau;
	double pairmu[3]={0};
	double pairen[3][3]={0};
	double mun;
	int i,j,k,l;
	
	/*declaration of NI variables*/
	int F=0;
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
	double mu,energy,expe,energyx;
	
	double p1a[3];
	double p1b[3]={0};
	double p2a[3][3],p2b[3][3];
	double p4b[3][3][3][3];
	double p4a[3][3][3][3]={0};
	double deltaP[90];
	double delta=1;
	double *point,*pair;
	double Grandpotential[ITERATIONMAX]={0};

	struct Configuration3ele conf;
	double inipoint[4][3]={{0},{0},{0},{0}};

	/*initial variables*/
	if (elemu!=0)
		tau=TAU*0.01/(1000*pow(10,4.2)*fabs(elemu));
	else
		tau=TAU;
	
	pairmu[0]=*(paireE+0)+elemu;
	pairmu[1]=*(paireE+1);
	pairen[0][1]=*(paireE+2);
	pairen[1][0]=*(paireE+2);
	pairen[0][2]=*(paireE+3);
	pairen[2][0]=*(paireE+3);
	pairen[1][2]=*(paireE+4);
	pairen[2][1]=*(paireE+4);
	mun=*(paireE+5);
	
	//For 3 symmetries, use input pairprobablilities inipro to generate point probabilities 
	//initial p2a
	if(memcmp("A3B",symmetry,3)==0)
	{	
		for(i=0;i<3;i++)
		{	
			//sum along index 1
			inipoint[i][0]=inipro[0][0]+inipro[0][1]+inipro[0][2]; 
			inipoint[i][1]=inipro[0][3]+inipro[0][4]+inipro[0][5];
			inipoint[i][2]=inipro[0][6]+inipro[0][7]+inipro[0][8];
		}
		
		//sum along index 0			
		inipoint[3][0]=inipro[0][0]+inipro[0][3]+inipro[0][6];
		inipoint[3][1]=inipro[0][1]+inipro[0][4]+inipro[0][7];
		inipoint[3][2]=inipro[0][2]+inipro[0][5]+inipro[0][8];
		
		//initial p2a
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				p2a[i][j]=(3*inipro[0][3*i+j]+3*inipro[1][3*i+j])/6.0;
	}
	else if (memcmp("AB",symmetry,2)==0)
	{
		
		for(i=0;i<2;i++)
		{
			//sum along index 1
			inipoint[i][0]=inipro[1][0]+inipro[1][1]+inipro[1][2]; 
			inipoint[i][1]=inipro[1][3]+inipro[1][4]+inipro[1][5];
			inipoint[i][2]=inipro[1][6]+inipro[1][7]+inipro[1][8];
		}
		for(i=2;i<4;i++)
		{
			//sum along index 0
			inipoint[i][0]=inipro[1][0]+inipro[1][3]+inipro[1][6]; 
			inipoint[i][1]=inipro[1][1]+inipro[1][4]+inipro[1][7];
			inipoint[i][2]=inipro[1][2]+inipro[1][5]+inipro[1][8];
		}
	 
	 	//initial p2a
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				p2a[i][j]=(inipro[0][3*i+j]+4*inipro[1][3*i+j]+inipro[2][3*i+j])/6.0;
	}
	else 
	{
		for(i=0;i<4;i++)
		{
			//sum along index 1
			inipoint[i][0]=inipro[0][0]+inipro[0][3]+inipro[0][6]; 
			inipoint[i][1]=inipro[0][1]+inipro[0][4]+inipro[0][7];
			inipoint[i][2]=inipro[0][2]+inipro[0][5]+inipro[0][8];
		}
		
		
		//initial p2a
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				p2a[i][j]=inipro[0][3*i+j];
	}
	
	
	/*natural iteration */
    while( delta > tau && F < ITERATIONMAX )
	{
		/*vaibles*/
		double lambda; 
		double p2x,p1x,f2x,f1x;
		double epoint=0;
		double epair=0;
		double ge=0;
		
		/*sum over exp(1/2 lambda beta)*/
		lambda=0;
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						mu=pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l];
						energy=0.5*(pairen[i][j]+pairen[i][k]+pairen[i][l]+pairen[j][k]+pairen[j][l]+pairen[k][l]);
						expe=2*beta*energy-1.5*beta*mu;
						
						if(expe>400)
							energyx=5e+173;
						else if(expe<-1e10)
							energyx=1.0+expe;
						else 
							energyx=exp(expe);							

						
						if(F==0)
						{
							if(memcmp("A3B",symmetry,3)==0)
							{
								p2x=pow(inipro[0][3*i+l]*inipro[0][3*j+l]*inipro[0][3*k+l]*inipro[1][3*i+j]*inipro[1][3*i+k]*inipro[1][3*j+k],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);									
							}
							else if(memcmp("AB",symmetry,2)==0)
							{
								p2x=pow(inipro[0][3*i+j]*inipro[1][3*i+k]*inipro[1][3*i+l]*inipro[1][3*j+k]*inipro[1][3*j+l]*inipro[2][3*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
							else 
							{
								p2x=pow(inipro[0][3*i+j]*inipro[0][3*i+k]*inipro[0][3*i+l]*inipro[0][3*j+k]*inipro[0][3*j+l]*inipro[0][3*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
						}
						else 
						{	
							f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
							f2x=p2aij[i][j]*p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l]*p2akl[k][l];
							p2x=pow(f2x,0.5);
							if(f1x<1e-70)
								p1x=5.6234132519e+43;
							else 
								p1x=pow(f1x,-0.625);
						}
							
						lambda+=energyx*p2x*p1x;
							
					}
				}
			}
		}
		/*iteration*/
		for(i=0;i<3;i++)
		{
			for(j=0;j<3;j++)
			{
				for(k=0;k<3;k++)
				{
					for(l=0;l<3;l++)
					{
						mu=pairmu[i]+pairmu[j]+pairmu[k]+pairmu[l];
						energy=0.5*(pairen[i][j]+pairen[i][k]+pairen[i][l]+pairen[j][k]+pairen[j][l]+pairen[k][l]);
						expe=2*beta*energy-1.5*beta*mu;
						
						if(expe>400)
							energyx=5e+173;
						else if(expe<-1e10)
							energyx=1.0+expe;
						else 
							energyx=exp(expe);							
				
						if(F==0)
						{
							if(memcmp("A3B",symmetry,3)==0)
							{
								p2x=pow(inipro[0][3*i+l]*inipro[0][3*j+l]*inipro[0][3*k+l]*inipro[1][3*i+j]*inipro[1][3*i+k]*inipro[1][3*j+k],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);									
							}
							else if(memcmp("AB",symmetry,2)==0)
							{
								p2x=pow(inipro[0][3*i+j]*inipro[1][3*i+k]*inipro[1][3*i+l]*inipro[1][3*j+k]*inipro[1][3*j+l]*inipro[2][3*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
							else 
							{
								p2x=pow(inipro[0][3*i+j]*inipro[0][3*i+k]*inipro[0][3*i+l]*inipro[0][3*j+k]*inipro[0][3*j+l]*inipro[0][3*k+l],0.5);
								f1x=inipoint[0][i]*inipoint[1][j]*inipoint[2][k]*inipoint[3][l];
								if(f1x<1e-70)
									p1x=5.6234132519e+43;
								else 
									p1x=pow(f1x,-0.625);								
							}
						}
						else 
						{	
							f1x=p1ai[i]*p1aj[j]*p1ak[k]*p1al[l];
							f2x=p2aij[i][j]*p2aik[i][k]*p2ail[i][l]*p2ajk[j][k]*p2ajl[j][l]*p2akl[k][l];
							p2x=pow(f2x,0.5);
							if(f1x<1e-70)
								p1x=5.6234132519e+43;
							else 
								p1x=pow(f1x,-0.625);
						}
							
						p4b[i][j][k][l]=energyx*p2x*p1x/lambda;
							
					}
				}
			}
		}
		
		/*generate 2 point site probablilities*/
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
			{	
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
		
		/*initial p2b*/
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				p2b[i][j]=(p2bij[i][j]+p2bik[i][j]+p2bil[i][j]+p2bjk[i][j]+p2bjl[i][j]+p2bkl[i][j])/6.0;
				
		
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
					for(l=0;l<3;l++)
						deltaP[9+3*(3*(3*i+j)+k)+l]=p4a[i][j][k][l]-p4b[i][j][k][l];
		
		/*caclulate the absnorm deltap*/
		delta=0;
		for(i=0;i<90;i++)
			delta+=fabs(deltaP[i]);
		
		
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
					for(l=0;l<3;l++)
						p4a[i][j][k][l]=p4b[i][j][k][l];
				
		/*grandpotential*/	
		for(i=0;i<3;i++)
			{
				epoint+=p1b[i]*pairmu[i];
			}
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				{
					epair+=0.5*p2b[i][j]*pairen[i][j];
				}
				/*entrop log of number of macro state*/
		ge=0;
		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
			{
				if(p2bij[i][j]>=LOGERR)
				{
					if(p2bik[i][j]>=LOGERR)
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bik[i][j]*(log(p2bik[i][j])-1));
								}
							}
						}
					}
					else //ik
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bij[i][j]*(log(p2bij[i][j])-1));
								}
							}
						}
					}
				}
				else // ij
				{
					if(p2bik[i][j]>=LOGERR)
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bik[i][j]*(log(p2bik[i][j])-1));
								}
							}
						}
					}
					else //ik
					{
						if(p2bil[i][j]>=LOGERR)
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else //kl
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else //jl
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else //kl
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bil[i][j]*(log(p2bil[i][j])-1));
								}
							}
						}
						else //il
						{
							if(p2bjk[i][j]>=LOGERR)
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjk[i][j]*(log(p2bjk[i][j])-1));
								}
							}
							else
							{
								if(p2bjl[i][j]>=LOGERR)
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bjl[i][j]*(log(p2bjl[i][j])-1)+p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=(p2bjl[i][j]*(log(p2bjl[i][j])-1));
								}
								else
								{
									if(p2bkl[i][j]>=LOGERR)
										ge+=(p2bkl[i][j]*(log(p2bkl[i][j])-1));
									else
										ge+=0;
								}
							}
						}
					}
				}	
			}
		for(i=0;i<3;i++)
		{
			if(p1bi[i]>=LOGERR)
			{
				if(p1bj[i]>=LOGERR)
				{
					if(p1bk[i]>=LOGERR)
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1));	
					}
					else
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bj[i]*(log(p1bj[i])-1));
					}
				}
				else
				{
					if(p1bk[i]>LOGERR)
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bk[i]*(log(p1bk[i])-1));
					}
					else
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bi[i]*(log(p1bi[i])-1));	
					}
				}
			}
			else
			{
				if(p1bj[i]>=LOGERR)
				{
					if(p1bk[i]>=LOGERR)
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bk[i]*(log(p1bk[i])-1));	
					}
					else
					{
						if(p1bl[i]>=LOGERR)
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bj[i]*(log(p1bj[i])-1));
					}
				}
				else
				{
					if(p1bk[i]>LOGERR)
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bk[i]*(log(p1bk[i])-1)+p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*(p1bk[i]*(log(p1bk[i])-1));
					}
					else
					{
						if(p1bl[i]>LOGERR)
							ge+=-1.25*(p1bl[i]*(log(p1bl[i])-1));
						else
							ge+=-1.25*0;	
					}
				}
			}
		}

		for(i=0;i<3;i++)
			for(j=0;j<3;j++)
				for(k=0;k<3;k++)
					for(l=0;l<3;l++)
					{	
						if(p4b[i][j][k][l]>LOGERR)
							ge+=-2*p4b[i][j][k][l]*(log(p4b[i][j][k][l])-1);			
						else 
							ge+=0;
					}
			
		
		Grandpotential[F-1]=epoint*Z-2*(epair)*Z-mun*Z-KB*tem*ge;
		
	}
	/*initial return value*/	
	memcpy(conf.Atoms,atoms,6);
	conf.NumIt=F;
	conf.Tem=tem;
	conf.Mu=elemu;
	conf.Grandenergy=Grandpotential[F];
	conf.pointP[0]=p1bi[0];
	conf.pointP[1]=p1bi[1];
	conf.pointP[2]=p1bi[2];
	conf.pointP[3]=p1bj[0];
	conf.pointP[4]=p1bj[1];
	conf.pointP[5]=p1bj[2];
	conf.pointP[6]=p1bk[0];
	conf.pointP[7]=p1bk[1];
	conf.pointP[8]=p1bk[2];
	conf.pointP[9]=p1bl[0];
	conf.pointP[10]=p1bl[1];
	conf.pointP[11]=p1bl[2];
	
	return conf;
}

void print3ele(struct Configuration3ele conf)
{
	int i;
	printf("Conf Atoms: %s\n",conf.Atoms);
	printf("iteration number is: %d\n",conf.NumIt);
	printf("the temperature is %.2f\n",conf.Tem);
	printf("the chemical potential is %.4f\n",conf.Mu);
	printf("grand potential is %.4f\n",conf.Grandenergy);
	printf("the probablility is:\n");
	for(i=0;i<4;i++)
		printf("%8.3f%8.3f%8.3f\n",conf.pointP[3*i],conf.pointP[3*i+1],conf.pointP[3*i+2]);
	
}

void phasedata2eledisorder(char *atoms,double *paireE, double *tem, double *mu)
{
	double index1,index2;
	double inipro[3][4]={{0},{0},{0}};
	double p[NoPROBCONDITION][4]={{0.3,0.15,0.15,0.4},{0.25,0.25,0.25,0.25},{0.7,0.1,0.1,0.1},{0.4,0.05,0.05,0.5},{0.1,0.15,0.15,0.6}};
	struct Configuration2ele conf;
	double energy[NoPROBCONDITION];
	int Number[NoPROBCONDITION];
	double point[NoPROBCONDITION][8];
	int i,j,k;
	int minindex;
	int phase;
	FILE *fp;

	fp=fopen("B:/phasediagram/disorder6.txt","w");
	for(index1=tem[0];index1<tem[1];index1+=tem[2])
 		for(index2=mu[0];index2<mu[1];index2+=mu[2])
 			{
 			//use different initial probabilities to NI calculate 
			 for(i=0;i<NoPROBCONDITION;i++)
 				{	
 					//initial initial probabilities
 					for(j=0;j<4;j++)
 					{
						inipro[0][j]=p[i][j]; 
					}
					
					conf=NI2ele(atoms,paireE,inipro,"disorder",index1,index2);

					//store the result 
					energy[i]=conf.Grandenergy;
					Number[i]=conf.NumIt;
					for(j=0;j<8;j++)
					{
						point[i][j]=conf.pointP[j];
					}

					
				}
			//give the minimum grandpotential index;
			minindex=0;
			for(i=1;i<NoPROBCONDITION;i++)
				{	
					if(energy[minindex]>=energy[i])
						minindex=i;
				}
			
			phase=phase2elelabel(point[minindex]);
			fprintf(fp,"%.0f %6.5f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d %10.6f %d\n",index1,index2,point[minindex][0],point[minindex][1],point[minindex][2],point[minindex][3],point[minindex][4],point[minindex][5],point[minindex][6],point[minindex][7],Number[minindex],energy[minindex],phase);
			}
			
	fclose(fp);
}

void phasedata2eleAB(char *atoms,double *paireE, double *tem,double *mu)
{	
	/*tem, mu  three numbers are min,max, decrease*/
	/*five different probabilities */
	double px[NoPROBCONDITION][4]={{0.4,0.1,0.1,0.4},{0.2,0.4,0.4,0.0},{0.2,0.3,0.3,0.2},{0.7,0.1,0.1,0.1},{0.8,0.05,0.05,0.1}};
	double py[NoPROBCONDITION][4]={{0.3,0.0,0.1,0.6},{0.3,0.6,0.1,0.0},{0.5,0.1,0.1,0.3},{0.25,0.25,0.25,0.25},{0.8,0.1,0.05,0.05}};
	double pz[NoPROBCONDITION][4]={{0.2,0.05,0.05,0.7},{0.2,0.3,0.3,0.2},{0.1,0.2,0.2,0.5},{0.6,0.15,0.15,0.1},{0.2,0.4,0.4,0.0}};
	double iniprobability[3][4]={{0},{0},{0}};
	int i,j,k;
	double index1,index2;
   
	struct Configuration2ele confAB;
	double energy[NoPROBCONDITION];
	double point[NoPROBCONDITION][8];
	int Num[NoPROBCONDITION];
	int minindex;
	int phase=0;
	FILE *fp;

	
	fp=fopen("B:/phasediagram/AB6.txt","w");
	for(index1=tem[0];index1<tem[1];index1+=tem[2])
		for(index2=mu[0];index2<mu[1];index2+=mu[2])
		{
			for(i=0;i<NoPROBCONDITION;i++)
			{	
		
				/*initial iniprobability*/
				for(j=0;j<4;j++)
				{
					iniprobability[0][j]=px[i][j];
					iniprobability[1][j]=py[i][j];
					iniprobability[2][j]=pz[i][j];
				}
	    		
				/*cacluate the ni*/
				confAB=NI2ele(atoms, paireE, iniprobability, "AB", index1, index2);				 
			 	/*store the iteration result*/				
				energy[i]=confAB.Grandenergy;
				Num[i]=confAB.NumIt;
				
				for(k=0;k<8;k++)
				{
					point[i][k]=confAB.pointP[k];
				}

				
			}
			
			/*grand energy's minimum index */
			minindex=0;
			for(i=1;i<NoPROBCONDITION;i++)
			{	
				if(energy[minindex]>=energy[i])
					minindex=i;
			}
			
			phase=phase2elelabel(point[minindex]);
			fprintf(fp,"%.0f %6.5f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d %10.6f %d\n",index1,index2,point[minindex][0],point[minindex][1],point[minindex][2],point[minindex][3],point[minindex][4],point[minindex][5],point[minindex][6],point[minindex][7],Num[minindex],energy[minindex],phase);
		}
	
	fclose(fp);	
}
	
void phasedata2eleA3B(char *atoms,double *paireE, double *tem,double *mu)
{	
	/*tem, mu  three numbers are min,max, decrease*/
	/*five different probabilities */
	double px[NoPROBCONDITION][4]={{0.3,0.0,0.2,0.5},{0.3,0.6,0.1,0.0},{0.5,0.1,0.1,0.3},{0.25,0.25,0.25,0.25},{0.8,0.1,0.05,0.05}};
	double py[NoPROBCONDITION][4]={{0.1,0.25,0.25,0.4},{0.1,0.25,0.25,0.4},{0.1,0.2,0.2,0.5},{0.1,0.4,0.4,0.1},{0.2,0.4,0.4,0.0}};
	double iniprobability[3][4]={{0},{0},{0}};
	int i,j,k;
	double index1,index2;
   
	struct Configuration2ele confA3B;
	//struct Configuration2ele confarrayAB[NoPROBCONDITION];
	double energy[NoPROBCONDITION];
	double point[NoPROBCONDITION][8];
	int Num[NoPROBCONDITION];
	int minindex;
	int phase;
	FILE *fp;

	
	fp=fopen("B:/phasediagram/A3B6.txt","w");
	for(index1=tem[0];index1<tem[1];index1+=tem[2])
		for(index2=mu[0];index2<mu[1];index2+=mu[2])
		{
			for(i=0;i<NoPROBCONDITION;i++) 
			{	
				/*initial iniprobability*/
				for(j=0;j<4;j++)
				{
					iniprobability[0][j]=px[i][j];
					iniprobability[1][j]=py[i][j];

				}
	    		
				/*cacluate the ni*/
				confA3B=NI2ele(atoms, paireE, iniprobability, "A3B", index1, index2);	
	
				/*store the iteration result*/
				energy[i]=confA3B.Grandenergy;
				Num[i]=confA3B.NumIt;
				for(j=0;j<8;j++)
				{
					point[i][j]=confA3B.pointP[j];
				}
				
			}
			
			/*grand energy's minimum index */
			minindex=0;
			for(i=1;i<NoPROBCONDITION;i++)
			{
				if(energy[minindex]>=energy[i])
					minindex=i;
			}
				
				
			/*write to files */
			phase=phase2elelabel(point[minindex]);
			fprintf(fp,"%.0f %6.5f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %d %10.6f %d\n",index1,index2,point[minindex][0],point[minindex][1],point[minindex][2],point[minindex][3],point[minindex][4],point[minindex][5],point[minindex][6],point[minindex][7],Num[minindex],energy[minindex],phase);			
		}
	fclose(fp);	
}

struct Diffenergy energydifference2ele(char *atoms, double *paireE, char *s_1, char *s_2, double tem,double mu)
{	//calculate energydifference of symmetry s1 and s2 phases. at given tem and chemical potential. 
	//return energydifference and the phaselabels 
	
	//initial probabilites. 
	double px_a[NoPROBCONDITION][4]={{0.3,0.15,0.15,0.4},{0.25,0.25,0.25,0.25},{0.7,0.1,0.1,0.1},{0.4,0.05,0.05,0.5},{0.1,0.15,0.15,0.6}};
	
	double px_b[NoPROBCONDITION][4]={{0.4,0.1,0.1,0.4},{0.2,0.4,0.4,0.0},{0.2,0.3,0.3,0.2},{0.7,0.1,0.1,0.1},{0.8,0.05,0.05,0.1}};
	double py_b[NoPROBCONDITION][4]={{0.3,0.0,0.1,0.6},{0.3,0.6,0.1,0.0},{0.5,0.1,0.1,0.3},{0.25,0.25,0.25,0.25},{0.8,0.1,0.05,0.05}};
	double pz_b[NoPROBCONDITION][4]={{0.2,0.05,0.05,0.7},{0.2,0.3,0.3,0.2},{0.1,0.2,0.2,0.5},{0.6,0.15,0.15,0.1},{0.2,0.4,0.4,0.0}};
	
	double px_c[NoPROBCONDITION][4]={{0.3,0.0,0.2,0.5},{0.3,0.6,0.1,0.0},{0.5,0.1,0.1,0.3},{0.25,0.25,0.25,0.25},{0.8,0.1,0.05,0.05}};
	double py_c[NoPROBCONDITION][4]={{0.1,0.25,0.25,0.4},{0.1,0.25,0.25,0.4},{0.1,0.2,0.2,0.5},{0.1,0.4,0.4,0.1},{0.2,0.4,0.4,0.0}};
	
	//variables 
	double inipro_1[3][4]={{0},{0},{0}};
	double inipro_2[3][4]={{0},{0},{0}};
	double energy_1[NoPROBCONDITION];
	double energy_2[NoPROBCONDITION];
	int minindex_1;
	int minindex_2;
	double point_1[NoPROBCONDITION][8];
	double point_2[NoPROBCONDITION][8];
	struct Configuration2ele conf_1;
	struct Configuration2ele conf_2;
	struct Diffenergy energydiff;
	int i,j,k,l;
	
	//two symmetries NI reslut and record the difference between two energies
	for(i=0;i<NoPROBCONDITION;i++)
	{	
		//judge the symmetries and initial pro
		if(memcmp("A3B",s_1,3)==0 && memcmp("AB",s_2,2)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_c[i][j];
				inipro_1[1][j]=py_c[i][j];
			
				inipro_2[0][j]=px_b[i][j];
				inipro_2[1][j]=py_b[i][j];
				inipro_2[2][j]=pz_b[i][j];
			}
		}
		else if(memcmp("AB",s_1,2)==0 && memcmp("A3B",s_2,3)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_b[i][j];
				inipro_1[1][j]=py_b[i][j];
				inipro_1[2][j]=pz_b[i][j];
				
				inipro_2[0][j]=px_c[i][j];
				inipro_2[1][j]=py_c[i][j];
			}
		}
		else if(memcmp("A3B",s_1,3)==0 && memcmp("dis",s_2,3)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_c[i][j];
				inipro_1[1][j]=py_c[i][j];
				
				inipro_2[0][j]=px_a[i][j];
			}
		}
		else if(memcmp("dis",s_1,3)==0 && memcmp("A3B",s_2,3)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_a[i][j];
				
				inipro_2[0][j]=px_c[i][j];
				inipro_2[1][j]=py_c[i][j];	
			}
		}
		else if(memcmp("dis",s_1,3)==0 && memcmp("AB",s_2,2)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_a[i][j];
					
				inipro_2[0][j]=px_b[i][j];
				inipro_2[1][j]=py_b[i][j];
				inipro_2[2][j]=pz_b[i][j];	
			}
		}
		else if(memcmp("AB",s_1,2)==0 && memcmp("dis",s_2,3)==0)
		{
			for(j=0;j<4;j++)
			{
				inipro_1[0][j]=px_b[i][j];
				inipro_1[1][j]=py_b[i][j];
				inipro_1[2][j]=pz_b[i][j];	
					
				inipro_2[0][j]=px_a[i][j];
			}
		}
				
		//NI calculation 
		conf_1=NI2ele(atoms,paireE,inipro_1,s_1,tem,mu);
		conf_2=NI2ele(atoms,paireE,inipro_2,s_2,tem,mu);
		//store the iteration result 
		energy_1[i]=conf_1.Grandenergy;
		energy_2[i]=conf_2.Grandenergy;
		
		for(j=0;j<8;j++)
		{
			point_1[i][j]=conf_1.pointP[j];
			point_2[i][j]=conf_2.pointP[j];
		}
	}
	//minimum energy 			
	minindex_1=0;
	minindex_2=0;
	for(i=0;i<NoPROBCONDITION;i++)
	{
		if(energy_1[minindex_1]>=energy_1[i])
			minindex_1=i;
		if(energy_2[minindex_2]>=energy_2[i])
			minindex_2=i;
	}
			
	//give the difference and label the phases 
	energydiff.G_energy=energy_1[minindex_1]-energy_2[minindex_2];
	energydiff.phase_1=phase2elelabel(point_1[minindex_1]);
	energydiff.phase_2=phase2elelabel(point_2[minindex_2]);
	for(i=0;i<8;i++)
	{
		energydiff.p_1[i]=point_1[minindex_1][i];
		energydiff.p_2[i]=point_2[minindex_2][i];
	}
	
	return energydiff;
}

void bisectionroot(char *atoms,double *paireE, FILE *fp_1, FILE *fp_2,FILE *fp_3,  char *s_1, char *s_2,double *tem,double *mu)
{	//find root of diffenergy use the method of bisection 
	//fp_1 A3B and dis phaseboundary
	//fp_2 AB and dis phaseboundary
	
	
	//vaiables 
	double index1,index2;  //temperature and chemical potential
	double left,right,mid;
	struct Diffenergy energydiff_1,energydiff_2,energydiff_3;


//	sprintf(fname,"boundary%s-%s1.txt",s_1,s_2);
	for(index1=tem[0];index1<=tem[1];index1+=tem[2])
	{	
		
		for(index2=mu[0];index2<mu[1];index2+=mu[2]) 
		{
			left=index2;
			right=index2+mu[2];
						
						
			//method of bisection
			while(right-left>DIFFDELTA)
			{	
				mid=(left+right)/2.0;
				energydiff_1=energydifference2ele(atoms,paireE,s_1,s_2,index1,left);
				energydiff_2=energydifference2ele(atoms,paireE,s_1,s_2,index1,mid);
				energydiff_3=energydifference2ele(atoms,paireE,s_1,s_2,index1,right);
				
				if((energydiff_1.G_energy<0 && energydiff_2.G_energy>0) ||(energydiff_1.G_energy>=0 && energydiff_2.G_energy<0))
					right=mid;
				else if ((energydiff_2.G_energy<0 && energydiff_3.G_energy>0) ||(energydiff_2.G_energy>=0 && energydiff_3.G_energy<0))
					left=mid;
				else				
					break;
			}	
		
			//write to file //

			//the format is tem,mu, phase_1,phase_2, energydiff
			
			if((energydiff_2.phase_1!=energydiff_2.phase_2) && (energydiff_2.G_energy<0.001)) //the phases one both side are different 
			{
			//	fprintf(fp,"%.0f %.4f  %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.6f\n",index1,mid,energydiff_2.phase_1,energydiff_2.p_1[0],energydiff_2.p_1[1],energydiff_2.p_1[2],energydiff_2.p_1[3],energydiff_2.p_1[4],energydiff_2.p_1[5],energydiff_2.p_1[6],energydiff_2.p_1[7],energydiff_2.phase_2,energydiff_2.p_2[0],energydiff_2.p_2[1],energydiff_2.p_2[2],energydiff_2.p_2[3],energydiff_2.p_2[4],energydiff_2.p_2[5],energydiff_2.p_2[6],energydiff_2.p_2[7],energydiff_2.G_energy);
				fprintf(fp,"%.0f %.4f  %d %d %.8f\n",index1,mid,energydiff_2.phase_1,energydiff_2.phase_2,energydiff_2.G_energy);
			}
		}
	}

}

void test(char *atoms,double *paireE,double tem, double mu)
{	//at give temperature and chemical potentilal, use different iniprobabilities to Nautral iteration caluclation.
	
	double p[6][4]={{0.1,0.2,0.2,0.5},{0.25,0.25,0.25,0.25},{0.7,0.1,0.1,0.1},{0.4,0.05,0.05,0.5},{0.3,0.1,0.1,0.5},{0.5,0.2,0.2,0.1}};
	double px[NoPROBCONDITION][4]={{0.25,0.25,0.25,0.25},{0.2,0.4,0.4,0.2},{0.7,0.1,0.1,0.1},{0.8,0.05,0.05,0.1},{0.3,0.0,0.0,0.7}};
	double py[NoPROBCONDITION][4]={{0.3,0.3,0.0,0.4},{0.3,0.6,0.1,0.0},{0.5,0.1,0.2,0.2},{0.05,0.8,0.1,0.05},{0.1,0.3,0.2,0.4}};
	double pz[NoPROBCONDITION][4]={{0.7,0.1,0.1,0.1},{0.1,0.2,0.2,0.5},{0.6,0.15,0.15,0.1},{0.2,0.4,0.4,0.0},{0.5,0.1,0.1,0.3}};
	struct Configuration2ele conf,conf1;
	double energy[NoPROBCONDITION];
	double point[NoPROBCONDITION][8];
	int Num[NoPROBCONDITION];
	FILE *f1,*f2;
	double inipro[3][4]={0};
	double inipro1[3][4]={0};
	int i,j;
	int index;
	
	
	f1=fopen("test1.txt","w");
	f2=fopen("test2.txt","w");
	for(i=0;i<NoPROBCONDITION;i++)
	{
		//initial iniprobabilities;
		for(j=0;j<4;j++)
			inipro[0][j]=p[i][j];
			
		//NI calculation
		conf=NI2ele(atoms,paireE,inipro,"disorder",tem,mu);

		//write to file;
		fprintf(f1,"%.3f %.3f %.3f %.3f %.0f %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d %10.6f\n",inipro[0][0],inipro[0][1],inipro[0][2],inipro[0][3],conf.Tem,conf.Mu,conf.pointP[0],conf.pointP[1],conf.pointP[2],conf.pointP[3],conf.pointP[4],conf.pointP[5],conf.pointP[6],conf.pointP[7],conf.NumIt,conf.Grandenergy);	
		
		for(j=0;j<4;j++)
				{
					inipro1[0][j]=px[i][j];
					inipro1[1][j]=py[i][j];
					inipro1[2][j]=pz[i][j];
				}
		conf1=NI2ele(atoms,paireE,inipro1,"AB",tem,mu);
		Num[i]=conf1.NumIt;
		energy[i]=conf1.Grandenergy;
		fprintf(f2,"%.3f %.3f %.3f %.3f %.0f %6.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %4d %10.6f\n",inipro1[0][0],inipro1[0][1],inipro1[0][2],inipro1[0][3],conf.Tem,conf.Mu,conf.pointP[0],conf.pointP[1],conf.pointP[2],conf.pointP[3],conf.pointP[4],conf.pointP[5],conf.pointP[6],conf.pointP[7],conf.NumIt,conf.Grandenergy);	
		
	}
	fclose(f1);
	fclose(f2);
	
	index=0;
	for(i=0;i<NoPROBCONDITION;i++)
	{
		if(energy[i]<=energy[index])
			index=i;
	}
	printf("the minimum index is %d\n", index);
	

}

int phase2elelabel(double *p)
{	//phaselabel phase_type==1, disorder phse

	//           phase_type==2 , A3B
	//			 phsae_type==3,  AB
	//           phase_type==4,  ABC 
	
	int phase_type;
	double p1bi[2],p1bj[2],p1bk[2],p1bl[2];
	int i;
	
	//transform the variables 
	for(i=0;i<2;i++)
	{
		p1bi[i]=*(p+0);
		p1bj[i]=*(p+2);
		p1bk[i]=*(p+4);	
		p1bl[i]=*(p+6);
	}	
	
	//label the phases 
	if((fabs(p1bi[0]-p1bj[0])<ERRP)&&(fabs(p1bi[1]-p1bj[1])<ERRP))
    {
       if((fabs(p1bj[0]-p1bk[0])<ERRP)&&(fabs(p1bj[1]-p1bk[1])<ERRP))
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=1;
           }
           else
           {
               phase_type=2;
           }
       }
       else
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=3;
           }
           else
           {
               if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
               {
                   phase_type=2;
               }
               else
               {
                   phase_type=4;
               }
           }
       }
    }
    else
    {
       if((fabs(p1bj[0]-p1bk[0])<ERRP)&&(fabs(p1bj[1]-p1bk[1])<ERRP))
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=2;
           }
           else
           {
               if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
               {
                   phase_type=3;
               }
               else
               {
                   phase_type=4;
               }
           }
       }
       else
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
              if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
              {
                  phase_type=2;
              }
              else
              {
                  phase_type=4;
              }
           }
           else
           {
              if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
              {
                  phase_type=4;
              }
              else
              {
                  if((fabs(p1bi[0]-p1bk[0])<ERRP)&&(fabs(p1bi[1]-p1bk[1])<ERRP))
                  {
                      if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
                      {
                          phase_type=3;
                      }
                      else
                      {
                          phase_type=4;
                      }
                  }
                  else
                  {
                      if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
                      {
                          phase_type=4;
                      }
                  }
              }
           }
       }
    }
    
	
	
	    
    return phase_type;
}

int phase3elelabel(double *p)
{	//phaselabel phase_type==1, disorder phse

	//           phase_type==2 , A3B
	//			 phsae_type==3,  AB
	//           phase_type==4,  ABC 
	
	int phase_type;
	double p1bi[3],p1bj[3],p1bk[3],p1bl[3];
	int i;
	
	for(i=0;i<3;i++)
	{
		p1bi[i]=*(p+0);
		p1bj[i]=*(p+3);
		p1bk[i]=*(p+6);	
		p1bl[i]=*(p+9);
	}
	
	
	if((fabs(p1bi[0]-p1bj[0])<ERRP)&&(fabs(p1bi[1]-p1bj[1])<ERRP))
    {
       if((fabs(p1bj[0]-p1bk[0])<ERRP)&&(fabs(p1bj[1]-p1bk[1])<ERRP))
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=1;
           }
           else
           {
               phase_type=2;
           }
       }
       else
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=3;
           }
           else
           {
               if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
               {
                   phase_type=2;
               }
               else
               {
                   phase_type=4;
               }
           }
       }
    }
    else
    {
       if((fabs(p1bj[0]-p1bk[0])<ERRP)&&(fabs(p1bj[1]-p1bk[1])<ERRP))
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
               phase_type=2;
           }
           else
           {
               if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
               {
                   phase_type=3;
               }
               else
               {
                   phase_type=4;
               }
           }
       }
       else
       {
           if((fabs(p1bk[0]-p1bl[0])<ERRP)&&(fabs(p1bk[1]-p1bl[1])<ERRP))
           {
              if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
              {
                  phase_type=2;
              }
              else
              {
                  phase_type=4;
              }
           }
           else
           {
              if((fabs(p1bi[0]-p1bl[0])<ERRP)&&(fabs(p1bi[1]-p1bl[1])<ERRP))
              {
                  phase_type=4;
              }
              else
              {
                  if((fabs(p1bi[0]-p1bk[0])<ERRP)&&(fabs(p1bi[1]-p1bk[1])<ERRP))
                  {
                      if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
                      {
                          phase_type=3;
                      }
                      else
                      {
                          phase_type=4;
                      }
                  }
                  else
                  {
                      if((fabs(p1bj[0]-p1bl[0])<ERRP)&&(fabs(p1bj[1]-p1bl[1])<ERRP))
                      {
                          phase_type=4;
                      }
                  }
              }
           }
       }
    }
    
	
	
	    
    return phase_type;
}







