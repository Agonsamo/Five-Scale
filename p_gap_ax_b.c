/*******************************************************************/
/*  p_gap_ax_b.c ->void P_GAP_AX_B(in,out,choice)                  */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the gap fraction in one crown       */
/*  Latest update                                 February 1997    */
/*******************************************************************/
/* inputs: view of solar zenith angle, LAI density, clumping and G */
/*******************************************************************/
                          
# include <math.h>
# include <stdio.h>
# include <string.h>
# include "data.h"
 
 
void P_GAP_AX_B(in_p,out_p,CHOICE)
 
struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;
 


{
double exp();
double ZA;
double Gs; 
out_p->Error[19] =0;
   
if(!strcmp("VZA",CHOICE))
{
	ZA=out_p->vza; 
	out_p->Gv=(in_p.A*ZA+in_p.C)*in_p.OMEGA_E/in_p.GAMMA_E;
	out_p->PgapV = exp(-out_p->Gv*in_p.LAI*in_p.B/(in_p.D*out_p->Vg_0*cos(ZA)));
	out_p->GFoliage = (in_p.A*ZA+in_p.C);
}
else if(!strcmp("SZA",CHOICE))
{
	ZA=in_p.SZA; 
	out_p->Gs=(in_p.A*ZA+in_p.C)*in_p.OMEGA_E/in_p.GAMMA_E;
	out_p->PgapS = exp(-out_p->Gs*in_p.LAI*in_p.B/(in_p.D*out_p->Sg_0*cos(ZA)));
	out_p->Lo = in_p.LAI*in_p.B/(in_p.D*out_p->Sg_0*cos(ZA));

}
else if(!strcmp("0",CHOICE))
{
	ZA=0; 
	Gs=(in_p.A*ZA+in_p.C)*in_p.OMEGA_E/in_p.GAMMA_E;
	out_p->Pgap0 = exp(-Gs*in_p.LAI*in_p.B/(in_p.D*in_p.R*in_p.R*PI*cos(ZA)));
}

else if(!strcmp("LAI",CHOICE))
{
	ZA=acos(0.537+0.025*in_p.LAI); 
	Gs=(in_p.A*ZA+in_p.C)*in_p.OMEGA_E/in_p.GAMMA_E;
	out_p->PgapV_mean = exp(-Gs*in_p.LAI*in_p.B/(in_p.D*out_p->Vg_0_mean*cos(ZA)));
}

}  
