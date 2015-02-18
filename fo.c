/*******************************************************************/
/*  fo.c ->void FO(in,out,choice)                                  */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  This function scale the gap fraction with vertical overlap     */
/*  or not.                                                        */
/*  Latest update                               October 27, 1999   */
/*******************************************************************/

 

# include <stdio.h>
# include <math.h>
# include <string.h>
# include "data.h"

void FO(in_p,out_p,CHOICE) 

struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;


{
	double exp();
	double SVG0=0,SVG=0;
	double PSVG=0,PSVG0=0;
	double theta=0;
 
 
	out_p->Fo=0;
	SVG0 = in_p.R*in_p.R*PI;
	PSVG0 = out_p->PSG_HOT0;
	if(!strcmp("VZA",CHOICE))
	{
		SVG   = out_p->Vg_0;
		PSVG  = out_p->PSG0_VIEW;
		theta = out_p->vza;
	}
	if(!strcmp("SZA",CHOICE))
	{
		SVG   = out_p->Sg_0;
		PSVG  = out_p->PSG0_SUN;
		theta = in_p.SZA;
	}
   
	if (SVG0 != 0) out_p->Fo = in_p.Fr*exp(-(SVG-SVG0)/SVG0*2*theta/PI); /* 08 October 1998 */

	out_p->Error[7] =0;

	if (out_p->Fo <0)
	{
		printf("\n Possible problem with Fo in sub FO()\n");
		out_p->Fo =0;
		out_p->Error[7] =1;
	}  
	if(!strcmp("VZA",CHOICE))
	{
		out_p->PVG = out_p->Pvg;
		out_p->Pvg = out_p->Pvg + (PSVG0-out_p->PVG_NADIR)*out_p->Fo; 
		out_p->Pv = out_p->Pv*(1-out_p->Fo);
	}
	if(!strcmp("SZA",CHOICE))
	{
		out_p->PIG = out_p->Pig;
		out_p->Pig = out_p->Pig + (PSVG0-out_p->PVG_NADIR)*out_p->Fo;
		out_p->Ps = out_p->Ps*(1-out_p->Fo);
	}


}  
