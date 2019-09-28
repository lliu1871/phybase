#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <R.h>

double caseIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1);
double caseIIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1);
double p0(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p1(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p2(double gamma0, double gamma1, double theta0, double theta1, int tree);
double p4(double gamma0, double gamma1, double theta0, double theta1, int tree);
double pxxx(double gamma0, double gamma1, double theta0, double theta1);
double pxxy(double gamma0, double gamma1, double theta0, double theta1);
double pyxx(double gamma0, double gamma1, double theta0, double theta1);
double pxyx(double gamma0, double gamma1, double theta0, double theta1);
double pxyz(double gamma0, double gamma1, double theta0, double theta1);
void tripleProb(double *gamma0, double *gamma1, double *theta0, double *theta1, double *p0,double *p1,double *p2,double *p3,double *p4);
double findMaxtrixElement(int row, int col, int nrow, double *a);
double rndu (void);
int nextNucleotide (int nucleotide, double branch);
void simnucleotides (int *nucleotide, int *length, double *branch, int *next);
void tripleLike1 (int *nuc, double *tau, double *theta, double *pi, double *likelihood) ;
void GetRNGstate(void);
void PutRNGstate(void);

int main(void)
{
	int i, nucleotide[1000], length=1000, next[1000];
	double branch=1.0;
	int nuc[5]={4,4,4,4,3};
	double tau[4]={0.01,0.001,0.003,0.002}, theta=0.001, pi[4]={0.25,0.25,0.25,0.25},likelihood;

    GetRNGstate();
	/*for(i=0; i<1000; i++)
	{
		nucleotide[i] = i % 4 + 1;
		printf("n %d\n",nucleotide[i]);
	}*/
	
	tripleLike1 (nuc, tau, &theta, pi, &likelihood);
	//printf("like %f\n",likelihood);
	//exit(-1);
	simnucleotides(nucleotide, &length, &branch, next);
    PutRNGstate();
    
	for(i=0; i<1000; i++)
	{
		//printf("n %d %d\n",nucleotide[i], next[i]);
	}
	return(1);
}

/*nucleotide: 0 A 1 G 2 C 3 T*/
void simnucleotides (int *nucleotide, int *length, double *branch, int *next)
{
	int i;

	for(i=0; i< (*length); i++)
	{
		next[i] = nextNucleotide(nucleotide[i], (*branch));
	}

}



int nextNucleotide (int nucleotide, double branch)
{
	double prob[4], ran = rndu();
	int i;


	for(i=0; i<4; i++)
	{
		if(i == (nucleotide-1))
			prob[i] = 0.25*(1+3*exp(-4*branch/3));
	 	else
			prob[i] = 0.25*(1-exp(-4*branch/3));
	}

	if(ran < prob[0])
		return(1);
	else if(ran < prob[0] + prob[1])
		return(2);
	else if(ran < prob[0] + prob[1] + prob[2])
		return(3);
	else 
		return(4);
		
}	

double findMaxtrixElement(int row, int col, int nrow, double *a)
{
	return (a[nrow*(col-1) + row-1]);
}

double rndu (void)
{
	double random = rand();
	return(random);
}

void tripleProb(double *gamma0, double *gamma1, double *theta0, double *theta1, double *p0,double *p1,double *p2,double *p3,double *p4)
{ 
	
	*p0 = pxxx(*gamma0,*gamma1,*theta0,*theta1);
	*p1 = pxxy(*gamma0,*gamma1,*theta0,*theta1);
	*p2 = pyxx(*gamma0,*gamma1,*theta0,*theta1);
	*p3 = *p2;
	*p4 = pxyz(*gamma0,*gamma1,*theta0,*theta1);

}



double caseIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1)
{
	double p;
    	p = 4*(1-exp((a1-2/theta1)*gamma0))/(theta0*theta1*(-a0+2/theta0)*(-a1+2/theta1));
    	return(p);
}

double caseIIprob(double gamma0, double gamma1, double theta0, double theta1, double a0, double a1)
{
	double p;
    	p = 12.0/(theta0*theta0*(-a0+2/theta0)*(-a1+6/theta0));
    	return(p);
} 




double p0(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = phy + 3*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8.0*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 6*exp(-(8.0*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 1 + 3*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8.0*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 6*exp(-(8.0*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
    	return (p);
}


double p1(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 3*phy + 9*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 6*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 3 + 9*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 6*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
//printf("P %f %f %f %f %d %f %f %f\n",gamma0,  gamma1, theta0,  theta1, tree, a0,  a1, p);
    	return (p);

}

double p2(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 3*phy - 3*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 3 - 3*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) + 6*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) - 6*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;

    	return (p);

}



double p4(double gamma0, double gamma1, double theta0, double theta1, int tree)
{
	double phy, a0, a1, p;

    	phy = 1-exp(-2*gamma0/theta1);
    	a0 = -8.0/3;
    	a1 = -8.0/3;

    	if(tree == 0)
    		p = 6*phy - 6*exp(-8.0/3*gamma1)*caseIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 12*exp(-8*(gamma0+gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 12*exp(-(8*gamma0+12*gamma1)/3)*caseIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
    	else
    	{
		p = 6 - 6*exp(-8.0/3*gamma1)*caseIIprob(gamma0,gamma1,theta0,theta1, 0, a1) - 12*exp(-8*(gamma0+gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, 0) + 12*exp(-(8*gamma0+12*gamma1)/3)*caseIIprob(gamma0,gamma1,theta0,theta1, a0, a1/2);
        	p = (1-phy)/3 * p;
    	}
    	p = p/16;
   
    	return (p);

}

double pxxx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
   	p = p0(gamma0, gamma1, theta0, theta1, 0) + 3 * p0(gamma0, gamma1, theta0, theta1, 1);
   	return (p);
}


double pxxy(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;   
	p = p1(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1); 
   	return(p);
}


double pyxx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;   
	p = p2(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1);  
   	return(p);
}


double pxyx(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
   	p = p2(gamma0, gamma1, theta0, theta1, 0) + p1(gamma0, gamma1, theta0, theta1, 1) + 2*p2(gamma0, gamma1, theta0, theta1, 1); 
   	return(p);
}


double pxyz(double gamma0, double gamma1, double theta0, double theta1)
{
	double p;
	p = p4(gamma0, gamma1, theta0, theta1, 0) + 3*p4(gamma0, gamma1, theta0, theta1, 1);
   	return(p);
}

void tripleLike1 (int *nuc, double *tau, double *theta, double *pi, double *likelihood) 
{

	int deta[4], index,i, j, k;
	double x1_a, x1_0, x1_1, x1_2, x1_3, x2_0, x2_1, x2_2, x2_3;
	double likelihood1, likelihood2, likelihood3, likelihood4;
	double f[10];

	//printf("%d %d %d %d %d %f %f %f %f %f %f %f %f %f\n",nuc[0],nuc[1],nuc[2],nuc[3],nuc[4],tau[0],tau[1],tau[2],tau[3],*theta,pi[0],pi[1],pi[2],pi[3]);
/*#########################################################
#
#                  caseI
#########################################################*/
deta[0] = 0;
deta[1] = 0;
deta[2] = 0;
deta[3] = 0;

if(nuc[1]!=nuc[2]) deta[0]=1;
if(nuc[1]!=nuc[3]) deta[1]=1;
if(nuc[0]!=nuc[4]) deta[2]=1;
if(nuc[0]!=nuc[1]) deta[3]=1;

x1_a = 3/(3-2*(*theta))*(1-exp(-(2/(*theta)-4.0/3)*tau[3]));
x1_0 = 1-exp(-2*tau[3]/(*theta));
x1_1 = 3/(3+2*(*theta))*(1-exp(-(2/(*theta)+4.0/3)*tau[3]));
x1_2 = 3/(3+4*(*theta))*(1-exp(-(2/(*theta)+8.0/3)*tau[3]));
x1_3 = 1/(1+2*(*theta))*(1-exp(-(2/(*theta)+4)*tau[3]));

x2_0 = 1;
x2_1 = 3/(3+2*(*theta));
x2_2 = 3/(3+4*(*theta));
x2_3 = 1/(1+2*(*theta));


likelihood1 = 0.00390625 * x1_0 * x2_0;

f[0] = x1_1 * x2_0 * exp(-4.0/3*(tau[0]));
f[1] = x1_1 * x2_0 * exp(-4.0/3*(tau[1]));
f[2] = x1_0 * x2_1 * exp(-4.0/3*(tau[2]));
f[3] = x1_a * x2_1 * exp(-4.0/3*(tau[3]));


for(i=0; i<4; i++)
{
	likelihood1 += 0.015625*(0.75-deta[i])*f[i];
}

f[0] = x1_2 * x2_0 * exp(-4.0/3*(tau[0]+tau[1]));
f[1] = x1_1 * x2_1 * exp(-4.0/3*(tau[0]+tau[2])); 
f[2] = x1_0 * x2_1 * exp(-4.0/3*(tau[0]+tau[3]));
f[3] = x1_1 * x2_1 * exp(-4.0/3*(tau[1]+tau[2])); 
f[4] = x1_0 * x2_1 * exp(-4.0/3*(tau[1]+tau[3]));  
f[5] = x1_a * x2_2 * exp(-4.0/3*(tau[2]+tau[3]));  
  
index = 0;
for(i=0; i<3; i++)
	for(j=(i+1); j<4; j++)
	{
		likelihood1 += 0.0625*(0.75-deta[i])*(0.75-deta[j])*f[index];
		index++;
	}

index=0;

f[0] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]));
f[1] = x1_1 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+tau[3]));
f[2] = x1_0 * x2_2 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));
f[3] = x1_0 * x2_2 * exp(-4.0/3*(tau[1]+tau[2]+tau[3]));

for(i=0;i<2;i++)
	for(j=(i+1);j<3;j++)
		for(k=(j+1);k<4;k++)
		{
			likelihood1 += (0.25)*(0.75-deta[i])*(0.75-deta[j])*(0.75-deta[k])*f[index];
			index++;
		}

f[0] = x1_1 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+tau[3]));
likelihood1 += (0.75-deta[0])*(0.75-deta[1])*(0.75-deta[2])*(0.75-deta[3])*f[0];


/*###############################################################################################
#
#                  caseII ((x,y),z)
#
#############################################################################################
*/
deta[0] = 0;
deta[1] = 0;
deta[2] = 0;
deta[3] = 0;

if(nuc[1]!=nuc[2]) deta[0]=1;
if(nuc[1]!=nuc[3]) deta[1]=1;
if(nuc[0]!=nuc[4]) deta[2]=1;
if(nuc[0]!=nuc[1]) deta[3]=1;

x1_0 = 1;
x1_1 = 3/(3+2*(*theta));
x1_2 = 3/(3+4*(*theta));
x1_3 = 1/(1+2*(*theta));

x2_0 = 1;
x2_1 = 3/(3+2*(*theta));
x2_2 = 3/(3+4*(*theta));
x2_3 = 1/(1+2*(*theta));


likelihood2 = 0.00390625 * x1_0 * x2_0;

f[0] = x1_1 * x2_0 * exp(-4.0/3*(tau[0]+tau[3]));
f[1] = x1_1 * x2_0 * exp(-4.0/3*(tau[1]+tau[3]));
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[2]));
f[3] = x1_0 * x2_1;


for(i=0;i<4;i++)
{
	likelihood2 += 0.015625*(0.75-deta[i])*f[i];
}

f[0] = x1_2 * x2_0 * exp(-4.0/3*(tau[0]+tau[1]+2*tau[3]));
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));  
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[0]+tau[3]));
f[3] = x1_2 * x2_1 * exp(-4.0/3*(tau[1]+tau[2]+tau[3]));  
f[4] = x1_1 * x2_1 * exp(-4.0/3*(tau[1]+tau[3]));  
f[5] = x1_1 * x2_2 * exp(-4.0/3*(tau[2]));  
  
index = 0;

for(i=0;i<3;i++)
	for(j=(i+1);j<4;j++)
	{
		likelihood2 += 0.0625*(0.75-deta[i])*(0.75-deta[j])*f[index];
		index++;
	}


index=0;

f[0] = x1_3 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+2*tau[3]));
f[2] = x1_2 * x2_2 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));
f[3] = x1_2 * x2_2 * exp(-4.0/3*(tau[1]+tau[2]+tau[3]));

for(i=0;i<2;i++)
	for(j=(i+1);j<3;j++)
		for(k=(j+1);k<4;k++)
		{
			likelihood2 += (0.25)*(0.75-deta[i])*(0.75-deta[j])*(0.75-deta[k])*f[index];
			index++;
		}

f[0] = x1_3 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
likelihood2 += (0.75-deta[0])*(0.75-deta[1])*(0.75-deta[2])*(0.75-deta[3])*f[0];


/*####################################################################
#
#         caseII: ((x,z),y)
#
###################################################################
*/
deta[0] = 0;
deta[1] = 0;
deta[2] = 0;
deta[3] = 0;

if(nuc[1]!=nuc[2]) deta[0]=1;
if(nuc[1]!=nuc[4]) deta[1]=1;
if(nuc[0]!=nuc[3]) deta[2]=1;
if(nuc[0]!=nuc[1]) deta[3]=1;

x1_0 = 1;
x1_1 = 3/(3+2*(*theta));
x1_2 = 3/(3+4*(*theta));
x1_3 = 1/(1+2*(*theta));

x2_0 = 1;
x2_1 = 3/(3+2*(*theta));
x2_2 = 3/(3+4*(*theta));
x2_3 = 1/(1+2*(*theta));


likelihood3 = 0.00390625 * x1_0 * x2_0;

f[0] = x1_1 * x2_0 * exp(-4.0/3*(tau[0]+tau[3]));
f[1] = x1_1 * x2_0 * exp(-4.0/3*(tau[2]));
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[1]+tau[3]));
f[3] = x1_0 * x2_1;


for(i=0;i<4;i++)
{
	likelihood3 += 0.015625*(0.75-deta[i])*f[i];
}

f[0] = x1_2 * x2_0 * exp(-4.0/3*(tau[0]+tau[2]+tau[3])); 
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+2*tau[3]));  
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[0]+tau[3]));
f[3] = x1_2 * x2_1 * exp(-4.0/3*(tau[1]+tau[2]+tau[3]));  
f[4] = x1_1 * x2_1 * exp(-4.0/3*(tau[2]));  
f[5] = x1_1 * x2_2 * exp(-4.0/3*(tau[1]+tau[3]));  
  

index=0;
for(i=0;i<3;i++)
	for(j=(i+1);j<4;j++)
	{
		likelihood3 += 0.0625*(0.75-deta[i])*(0.75-deta[j])*f[index];
		index++;
	}


index=0;

f[0] = x1_3 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));
f[2] = x1_2 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+2*tau[3]));
f[3] = x1_2 * x2_2 * exp(-4.0/3*(tau[1]+tau[2]+tau[3]));

for(i=0;i<2;i++)
	for(j=(i+1);j<3;j++)
		for(k=(j+1);k<4;k++)
		{
			likelihood3 += (0.25)*(0.75-deta[i])*(0.75-deta[j])*(0.75-deta[k])*f[index];
			index++;
		}

f[0] = x1_3 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
likelihood3 += (0.75-deta[0])*(0.75-deta[1])*(0.75-deta[2])*(0.75-deta[3])*f[0];

/*###############################################################################################
#
#                  caseII ((z,y),x)
#
#############################################################################################
*/
deta[0] = 0;
deta[1] = 0;
deta[2] = 0;
deta[3] = 0;

if(nuc[1]!=nuc[4]) deta[0]=1;
if(nuc[1]!=nuc[3]) deta[1]=1;
if(nuc[0]!=nuc[2]) deta[2]=1;
if(nuc[0]!=nuc[1]) deta[3]=1;

x1_0 = 1;
x1_1 = 3/(3+2*(*theta));
x1_2 = 3/(3+4*(*theta));
x1_3 = 1/(1+2*(*theta));

x2_0 = 1;
x2_1 = 3/(3+2*(*theta));
x2_2 = 3/(3+4*(*theta));
x2_3 = 1/(1+2*(*theta));


likelihood4 = 0.00390625 * x1_0 * x2_0;

f[0] = x1_1 * x2_0 * exp(-4.0/3*(tau[2]));
f[1] = x1_1 * x2_0 * exp(-4.0/3*(tau[1]+tau[3]));
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[0]+tau[3]));
f[3] = x1_0 * x2_1;


for(i=0;i<4;i++)
{
	likelihood4 += 0.015625*(0.75-deta[i])*f[i];
}


f[0] = x1_2 * x2_0 * exp(-4.0/3*(tau[2]+tau[1]+tau[3]));
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));  
f[2] = x1_1 * x2_1 * exp(-4.0/3*(tau[2]));
f[3] = x1_2 * x2_1 * exp(-4.0/3*(tau[1]+tau[0]+2*tau[3]));  
f[4] = x1_1 * x2_1 * exp(-4.0/3*(tau[1]+tau[3]));  
f[5] = x1_1 * x2_2 * exp(-4.0/3*(tau[0]+tau[3]));  
  

index = 0;
for(i=0;i<3;i++)
	for(j=(i+1);j<4;j++)
	{
		likelihood4 += 0.0625*(0.75-deta[i])*(0.75-deta[j])*f[index];
		index++;
	}


index=0;

f[0] = x1_3 * x2_1 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
f[1] = x1_2 * x2_1 * exp(-4.0/3*(tau[2]+tau[1]+tau[3]));
f[2] = x1_2 * x2_2 * exp(-4.0/3*(tau[0]+tau[2]+tau[3]));
f[3] = x1_2 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+2*tau[3]));

for(i=0;i<2;i++)
	for(j=(i+1);j<3;j++)
		for(k=(j+1);k<4;k++)
		{
			likelihood4 += (0.25)*(0.75-deta[i])*(0.75-deta[j])*(0.75-deta[k])*f[index];
			index++;
		}

f[0] = x1_3 * x2_2 * exp(-4.0/3*(tau[0]+tau[1]+tau[2]+2*tau[3]));
likelihood4 = likelihood4 + (0.75-deta[0])*(0.75-deta[1])*(0.75-deta[2])*(0.75-deta[3])*f[0];

*likelihood = pi[nuc[0]-1] * (likelihood1 + exp(-2*tau[3]/(*theta))/3*(likelihood2+likelihood3+likelihood4));

//printf("likke %1.10f %1.10f %1.10f %1.10f %1.10f\n", likelihood1, likelihood2,likelihood3,likelihood4, *likelihood);
}

