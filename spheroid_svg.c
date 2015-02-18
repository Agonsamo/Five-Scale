/*******************************************************************/
/*  spheroid_svg.c ->void SPHEROID_SVG(in,out,choice)              */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates crown projection on the ground      */
/*  for deciduous (spheroid) shape                                 */    
/*  Latest update                                February 7, 1997  */
/*******************************************************************/


# include <math.h>
# include <stdio.h>
# include <string.h>
# include "data.h"

 
void SPHEROID_SVG(in_p,out_p,CHOICE)
 
struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;



{  

	double ZA;                             /* view or solar angle */

	double b=in_p.Hb/2.;              /* vertical axis of the spheroid */    
	double za_prime;      
	if(!strcmp("VZA",CHOICE)) ZA = out_p->vza;
	else if(!strcmp("SZA",CHOICE)) ZA = in_p.SZA;
	else if(!strcmp("LAI",CHOICE)) ZA =acos(0.537+0.025*in_p.LAI);

	if (ZA == PI/2.) ZA = ZA - 0.0000000000000000001; 
							/* to get an aproximated value when ZA = 90 deg */ 


	za_prime = atan(b/in_p.R*tan(ZA));

	if(!strcmp("VZA",CHOICE)) out_p->Vg_0 = PI*in_p.R*in_p.R/(cos(za_prime));
	else if(!strcmp("SZA",CHOICE)) out_p->Sg_0 = PI*in_p.R*in_p.R/(cos(za_prime));
	else if(!strcmp("LAI",CHOICE)) out_p->Vg_0_mean = PI*in_p.R*in_p.R/(cos(za_prime));




  
}
