/*******************************************************************/
/*  spheroid_ta.c ->void SPHEROID_TA(in,out)                       */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates crown projection to the viewer      */
/*  for deciduous (spheroid) shape                                 */ 
/*  Latest update                                February 1, 1997  */
/*******************************************************************/


# include <stdio.h>
# include <math.h>
# include <string.h>
# include "data.h"



void SPHEROID_TA(in_p,out_p)
 
struct PARAMETER in_p;
struct RESULT *out_p;


{  
  
	out_p->tab = PI*in_p.R*(in_p.Hb/2.*sin(out_p->vza)+in_p.R*cos(out_p->vza));   
	out_p->V = 2/3.*PI*in_p.R*in_p.R*in_p.Hb;  
    
}
      
