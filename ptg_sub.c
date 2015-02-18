/*******************************************************************/
/*  ptg_sub.c ->void PTG_SUB(in,out,choice,max,inc)                */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the hotspot kernels                 */
/*  Latest update                              September 25, 1998  */
/*******************************************************************/


# include <stdio.h>
# include <math.h>
# include <string.h>
# include "data.h"
 
void PTG_SUB(in_p,out_p,CHOICE,MAX,INC)
 
struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;
double INC;              /* increment of the numerical integration */
int MAX;                    /* maximum number of integration steps */ 

{
	double IN1 = 0;                   /* numerator value of integral */
	double IN2 = 0;                 /* denominator value of integral */
	double FLAG_IN1=0,FLAG_IN2=0;
	register int i=0;
	double exp();
	double F_theta=0;
	double Lt=0;
	double Cold=0,Hot=0;
	double W=0,PTG=0,F=0,H=0;
	double XI =0;
	double lambda_m=0;

	out_p->Error[10] =0;
	/**************calculation of F() ************************************/
	if(!strcmp("CANOPY",CHOICE))
	{

		W = in_p.Ws;
		Lt=out_p->Ls;
		Cold=out_p->PT_Cold;
		Hot = 1-out_p->Pig; 
		XI = out_p->xi;
		H = out_p->H;
		lambda_m = out_p->lambda_m;


	} else 
	if(!strcmp("GROUND",CHOICE))
	{
		W = out_p->Wt; 
		Lt=out_p->Lt;
		Cold=out_p->PS;
		Hot = out_p->Pig; 
		XI = out_p->xi;
		H = out_p->H;
		lambda_m = out_p->lambda_m;
	} else
	if(!strcmp("NADIR",CHOICE))
	{
		
		W = sqrt(PI*in_p.R*in_p.R);
		Lt = out_p->OmegaT*PI*in_p.R*in_p.R*in_p.D/in_p.B;
		Cold = out_p->Pvg*(1-out_p->Pig); 
		Hot = out_p->Viewed_shadow;
		out_p->Medium =   Hot - Cold;
		if(Hot - Cold < 0 )
		{
			Hot = Cold;
		}
		XI = out_p->vza;
		H=out_p->Hc/3. + in_p.Hb + in_p.Ha;         
		lambda_m =H*tan(XI); 
	}

	

	if (XI < PI/2.)
	{
		for(i=0;i<MAX ;i++)  
		{
			IN1 = IN1 + exp(-Lt*(1+(lambda_m + INC*i)/W))/atan((lambda_m +INC*i )/H) ;

			IN2 =IN2 + exp(-Lt*(1+(lambda_m + INC*i)/W));     


			if (FLAG_IN1 == IN1 &&FLAG_IN2 == IN2 ) i = MAX;
			FLAG_IN1 = IN1;
			FLAG_IN2 = IN2;
		} 


	if (IN1 > 10000 || IN1 < -10000)
		{   
			IN1 =0;
			if(XI>0.00001) 
			{
			/*	printf("\n (IN1) Possible problem with ");
				printf("Hotspot kernel at xi= %5.3f",
					XI*180./PI); */
				out_p->Error[10] =1;
			} 
		}

	if (IN2 > 10000 || IN2 < -10000) 
		{   
			IN2 =0;
			if(XI>0.00001)
			{
			/*	printf("\n (IN2) Possible problem with ");
				printf("Hotspot kernel at xi= %5.3f",XI*180/PI);*/
				out_p->Error[10] =1;
			}
	}   


	F_theta =0;

	if(IN2 !=0) F_theta = 1-IN1/IN2*XI; 


	F=F_theta;

	PTG=  (Hot-Cold)*F_theta + Cold; 

	}
	else 
	{  
		PTG = Cold;
		F =0;
	}  

       
       
        
	if (F_theta <0) PTG = Cold;
	if (F_theta >1) PTG = Cold;
	/***********************************************************************/

	if(!strcmp("CANOPY",CHOICE))
	{
		out_p->PT=PTG;
		out_p->Fs=F;
	} else
	if(!strcmp("GROUND",CHOICE))
	{
		out_p->PG=PTG;
		out_p->Ft=F;

	} else
	if(!strcmp("NADIR",CHOICE))
	{
		out_p->Error[12] =0;
	/*	if((out_p->Pvg> PTG) && (out_p->Viewed_shadow > out_p->Pvg-out_p->PS)) out_p->PS = (out_p->Pvg- PTG); */
		out_p->ZG = PTG;
		if(out_p->PS > out_p->Pvg - out_p->ZG)	out_p->PS = out_p->Pvg - out_p->ZG;
		if(out_p->Pvg - out_p->ZG < 0)
		{
			out_p->PS = 0;
			out_p->Error[12] = 1;
		}
		out_p->Fn = F*out_p->Fd;
	}
 
        

}
