/*******************************************************************/
/*  q.c ->void Q(in,out,DEL)                                       */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates Q1tot and Q2tot ...                 */    
/*  Latest update                                October 02, 2001  */
/*******************************************************************/

 
# include <math.h>
# include <stdio.h>
# include "data.h"


void Q(in_p,out_p,CP,shape)
 
struct PARAMETER in_p;
struct RESULT *out_p;
double CP;
int shape; // not used anymore


{
	double Q1,Q2;
	int i,j; 
	double a;
	double Pav_i;
	double DEL=0, DELBroadleaf=0;

	DEL = 1. - out_p->xi*CP/PI;
	if(DEL<0) DEL =0;
	if(DEL>1) DEL =1;


	DELBroadleaf = out_p->xi*CP/PI-CP*60/180;
	if(DELBroadleaf<0) DELBroadleaf =0;
	if(DELBroadleaf>1) DELBroadleaf =1;


    out_p->QQ1=0;

	out_p->QQ2=0; 

	out_p->Error[8] =0;


	Q1=(1.-exp(-(out_p->Cs*out_p->Lo_90+out_p->Cv*out_p->Lo_90))) *out_p->Cs*out_p->Cv/(out_p->Cs+out_p->Cv);
	if (Q1>1.) 
	{
		Q1=1.;
		out_p->Error[8] =1;
	}
	if (sqrt ((out_p->Cs - out_p->Cv)*(out_p->Cs - out_p->Cv))<0.0000001) out_p->Cs = out_p->Cs - 0.0000001; 

	Q2=(exp(-out_p->Cs*out_p->Lo_90)-exp(-out_p->Cv*out_p->Lo_90))*out_p->Cs*out_p->Cv/(out_p->Cv-out_p->Cs);
	if (Q2>1.) 
	{
		Q2=1.;
		out_p->Error[9] =1;	
	}

	a = out_p->Gv*in_p.LAI*in_p.B/(in_p.D*out_p->Vg_0*cos(out_p->vza))*cos(out_p->vza)/cos(in_p.SZA);

    for (i=1;i<NN;i++)
    { 
		Pav_i=0; 

		if (exp(-(i-1)*a) <0.00000000000000000000000001) i = NN; 
		for (j=i;j<NN;j++)
		{
			Pav_i = Pav_i + out_p->Ptreev[j] ;
			if(out_p->Ptreev[j] <0.00000000000000000000000001 && j>NN/2.) j=NN; /* allows very small probabilities ...*/
		}

		/* Pav_i represents the probability of going through the ith tree outline */ 
		out_p->QQ1 = out_p->QQ1 + 1.*Q1*Pav_i*pow(out_p->PgapV,i-1.)*exp(-(i-1)*a); /* times probability of having foliage at that height maybe */
																					/* PgapV should follow this too ... */
		out_p->QQ2 = out_p->QQ2 + 1.*Q2*Pav_i*pow(out_p->PgapV,i-1.)*exp(-(i-1)*a);

	} 


	if(out_p->QQ1>1)
	{
		out_p->QQ1=1.;
		out_p->Error[8] =1;
	}
	
	if(out_p->QQ2>1)
	{
		out_p->QQ2=1.;
		out_p->Error[9] =1;
	}

	if(in_p.GAMMA_E>1)
	{
		out_p->QQ1B = out_p->QQ1*(1-DEL);
		out_p->QQ2B = out_p->QQ2*(1-DEL);
	} else
	{
		out_p->QQ1B = out_p->QQ1*(DELBroadleaf);
		out_p->QQ2B = out_p->QQ2*(DELBroadleaf);
	}
	
 

    out_p->QQ1= out_p->QQ1*(DEL);
    out_p->QQ2= out_p->QQ2*(DEL);


}
 
 
 

