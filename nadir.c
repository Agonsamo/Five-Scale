/*******************************************************************/
/*  nadir.c ->void NADIR(in,out)                                   */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  This soubroutine works like the hotspot function, but is       */
/*  applied to the shaded background viewed. The original Pig*Pvg  */
/*  equation can induce and underestimation of shaded background   */
/*  Latest update                                September 1998    */
/*******************************************************************/

# include <math.h>
# include <stdio.h>
# include "data.h"
 
 
void NADIR(in_p,out_p)

struct PARAMETER in_p;
struct RESULT *out_p;

{
	double  Fd1 =0,Fd2=0, Fd3=0;
	double i=0;
	double Lt = 0;
	double Wt = 0;
	double dist=0;
	double total =0;

	/* recalculation of Lt and Wt */
	Lt = out_p->OmegaT*PI*in_p.R*in_p.R*in_p.D/in_p.B;
	Wt = PI*in_p.R*in_p.R;

	for(i=0.;i<100.;i=i+0.1 )
	{	
		dist = i*exp(-Lt*(1+i/Wt));
		if(i > out_p->H/tan(in_p.SZA)) 
		{
			Fd2 = Fd2 +exp(-Lt*(1+i/Wt));
			total = total + exp(-Lt*(1+i/Wt));
		}
		else 
		{
			Fd2 = Fd2+ dist/(out_p->H/tan(in_p.SZA));
		
			total =total + exp(-Lt*(1+i/Wt));
		}
	}

	Fd2 = Fd2/total; 

	if(in_p.Ha/tan(in_p.SZA) > 2*in_p.R  )  Fd1 =1;
	else Fd1 = in_p.Ha/(tan(in_p.SZA)*2*in_p.R);
	if(in_p.Ha/tan(in_p.SZA) > out_p->E_r) Fd1 =0;

	out_p->Error[11]=0;   
	Fd3= (PI*in_p.R*in_p.R*(Fd1-Fd2)+out_p->Sg_0*Fd2)/(out_p->Sg_0);
	if(Fd3 >1 ) 
	{
		Fd3 =1;
		out_p->Error[11]=1;
	}
	if(Fd3 <0)
	{
		Fd3 =0;
		out_p->Error[11]=1;
	}
	out_p->Viewed_shadow = (1-out_p->Pig)*Fd3 ; 	
	out_p->Fd = Fd3;
}

