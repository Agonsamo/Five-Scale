/*******************************************************************/
/*  cone_gs.c ->void CONE_GS(in,out,choice)                        */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates crown projection on the ground      */
/*  Latest update                                March 27, 1995    */
/*******************************************************************/
/*   3 possible cases:                                             */
/*   1 --> nadir (vza =0)                                          */
/*   2 --> half apex angle > vza > 0                               */
/*   3 --> vza >= half apex angle                                  */
/*******************************************************************/


# include <math.h>
# include <string.h>
# include "data.h"
 
void CONE_GS(in_p,out_p,CHOICE)
 
struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;


{  
	double gamma=0.;
	double g_sc=0.;
	double g_s=0.;
	double ZA; 

	if(!strcmp("VZA",CHOICE)) ZA=out_p->vza;
	else if(!strcmp("SZA",CHOICE)) ZA=in_p.SZA;
	else if(!strcmp("LAI",CHOICE)) ZA=acos(0.537+0.025*in_p.LAI);
    
    
   

	if (ZA==PI/2.)ZA=PI/2. - 0.000000000001; /* to get a value at za=90 deg */ 

   
	if (ZA == 0)                                                  /* 1 */
    {         
		g_s = PI*in_p.R*in_p.R;   
    }
	else if (ZA > 0 && ZA < in_p.ALPHA)                           /* 2 */
    {  
		g_s = 2*tan(ZA)*in_p.R*in_p.Hb +PI*in_p.R*in_p.R;
    }
	else if ((ZA >= in_p.ALPHA) && (ZA < PI/2))                   /* 3 */
	{  
		gamma= asin(tan(in_p.ALPHA)/tan(ZA));
		g_s  = tan(ZA)*in_p.Hb*in_p.R*2 + 
               (1./tan(gamma)+PI/2 +gamma)*in_p.R*in_p.R;
		g_sc = (1./tan(gamma)-PI/2 +gamma)*in_p.R*in_p.R;
	}
	if(!strcmp("VZA",CHOICE)) 
	{
		out_p->Vg_0=g_s;
		out_p->Vgc=g_sc;
	}
	else if(!strcmp("SZA",CHOICE)) 
	{
		out_p->Sg_0=g_s;
		out_p->Sgc=g_sc;
	}
	else if(!strcmp("LAI",CHOICE)) 
	{
		out_p->Vg_0_mean=g_s;
		out_p->Vgc=g_sc;
	}

  }

