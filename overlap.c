/*******************************************************************/
/*  overlap.c ->void OVERLAP(in,out,choice)                        */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the gap fraction (Pvg and Pig)      */
/*  and Pv or Ps  and Pvc or Psc                                   */
/*  Latest update                                 February 1997    */
/*******************************************************************/


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include "data.h"

void OVERLAP(in_p,out_p,CHOICE) 

struct PARAMETER in_p;  
struct RESULT *out_p; 
char *CHOICE;
{
	double PSV;                          /* overlap probability    */
	double PSVC;                              /* cone overlap prob */
	double Pgap;                         /* gap in one tree prob   */
	double ac;                                  /* cone projection */
	double A;                                      /* quadrat size */
	double Ptree[NN]; 
	double Ptreec[NN]; 
	double a;   /* crown projection on the ground, either Vg or Sg */
	double i_0;                /* mean number of trees per quadrat */ 


	double pow();
	int  v; 
	register int i,j; 
	double Ptj[NN],Pt;  
	double P;      
	double ii,jj;                
	double aa,PNN,Ptj1=0;
	double Ptjc[NN],Ptc,PC,Ptjc1=0; 
	double aac,PNNC;
 
	A=in_p.B/in_p.n; 
	i_0 = in_p.D/in_p.n;
	if(!strcmp("VZA",CHOICE))
	{
		Pgap = out_p->PgapV;
		a = out_p->Vg_0;
		ac= out_p->Vgc;  
	}
	else if(!strcmp("SZA",CHOICE))
	{
		Pgap = out_p->PgapS;
		a = out_p->Sg_0;
		ac= out_p->Sgc;  
	} else if (!strcmp("LAI",CHOICE)) /* used in multiple scattering view factors */
	{
		Pgap = out_p->PgapV_mean;
		a = out_p->Vg_0_mean;
		ac= 0;
	}
	else if (!strcmp("NADIR",CHOICE))
	{
		Pgap = out_p->Pgap0;
		a = PI*in_p.R*in_p.R;
		ac= 0;
	}

	P=0;
	PC=0;
	PSV =0;
	PSVC =0;


	for(j=0;j<NN;j++)
	{    
		Pt = 0;      /* initialisation of variables */ 
		Ptc =0;
		Ptree[j] =0;  
		Ptreec[j] =0;  
		for(i=j;i<NN;i++)
		{
			Ptj[i] = 0.;
			Ptjc[i] =0.;
			ii=i*1.;
			jj=j*1.;
			if(Ptj[i-1] <0.0000001 && ii>i_0)
			{
				Ptj[i]=0;
				Ptjc[i]=0;
			}
			/* missing an else if compared to Windows version */
			 else
			{
				aa =a; 
				aac=ac;    
				if (ii> i_0) aa= a*i_0/ii ; 
				if (a > A  ) aa =A;  
				if (ac > A  ) aac =A;  
				if (jj==0)
				{
					Ptj[i] = 0;
					if(ii!=0)Ptj[i]=out_p->Px[i]*pow(1.-aa/A,ii);
					Ptjc[i] = 0;
					if(ii!=0)Ptjc[i]=out_p->Px[i]*pow(1.-aac/A,ii);
				} else  
				{
					PNN = out_p->Px[i];
					PNNC =out_p->Px[i];
					if (aa !=A)  PNN  = out_p->Px[i]*pow(1.-aa/A,ii-jj); 
					if (aac !=A) PNNC = out_p->Px[i]*pow(1.-aac/A,ii-jj);
					for(v=1;v<=j;v++)
					{
						PNN = 1.*PNN*(i+j-v)*aa/(A*(v))*(1.-aa/A) ;
						PNNC = 1.*PNNC*(i+j-v)*aac/(A*(v))*(1.-aac/A) ;
					}
					Ptj[i] =PNN; 
					Ptjc[i]=PNNC;
               
				 }    
			}

			Ptj[0] = out_p->Px[0];           
			Ptree[0] = Ptj[0];         
			Pt = Ptj[i] +Pt;
			Ptjc[0] = out_p->Px[0];         
			Ptreec[0] = Ptjc[0];
			Ptc = Ptjc[i] +Ptc;
		} 
		if(j==1) 
		{ 
			Ptj1=Pt;
			Ptjc1=Ptc;
		}
		if((Pgap == 0) && (j == 0))
		{
			P = Pt;
			PC =Ptc;
		}else
		{
			Ptree[j] = Pt;                /* November 29, 1996 */  
			Ptreec[j] =Ptc; 
			/* for( )Pow */
			P = P + Pt*pow(Pgap,j*1.); 		/* it should be possible to speed up the process by checking when */
			PC = PC + Ptc*pow(Pgap,j*1.); 	/* Pt and/or pow(Pgap,j) goes near zero	 ... */

			/* if I replace Pgap with Pgap[j], it will add a height dependency */

			if(!strcmp("VZA",CHOICE))   
  			{
   				out_p->Ptreev[j]=Ptree[j];
			}
			if(!strcmp("SZA",CHOICE))
  			{
   				out_p->Ptrees[j]=Ptree[j];
			}
		}
        
		PSV = 1 - P -  Ptj1*(1-Pgap); 
		PSVC = 1 - PC -  Ptjc1*(1-Pgap);       
		PSV = (PSV*a-PSVC*ac)/(a-ac);        
	}

	/* here's a lot of probably useless check on the probability, except at
	   large vza  */

	if(P <=0 || P > 1)
	{
		P = 0;
		PSV = 1;
	}
	if(PC <=0 || PC > 1)
	{
		PC = 0;
		PSVC = 1;
	}
	if(a/A >=  1)
	{
		PSV = 1;
		P=0;
	}
	if(ac/A >=  1)
	{
		 PSVC = 1;
		 PC=0;
	}

	if (PSV >=1) PSV =1 ;

	if (P<0.00001) PSV = 1;
 
	if (PSVC >=1) PSVC =1 ;

	if (PC<0.00001) PSVC = 1;
	if (ac <= 0) PSVC =0;

	if(!strcmp("VZA",CHOICE))
	{
		out_p->Pvg=P;
		out_p->Pv=PSV;
		out_p->Pvc=PSVC;
	} else
	if(!strcmp("SZA",CHOICE))
	{
		out_p->Pig=P;
		out_p->Ps=PSV;
		out_p->Psc=PSVC;
	} else
	if(!strcmp("LAI",CHOICE))
	{
		out_p->Pvg_mean=P;
	}
	else
	if(!strcmp("NADIR",CHOICE))
	{
		out_p->PVG_NADIR=P;
	}


} 
