/*******************************************************************/
/*  distance.c ->void DISTANCE(in,out)                             */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the mean gap between tree crowns    */
/*  Latest update                              November, 1999      */
/*******************************************************************/


# include <stdio.h>
# include <math.h>
# include <string.h>
# include "data.h"

void DISTANCE(in_p,out_p) 

struct PARAMETER in_p;
struct RESULT *out_p;

{
  
	double Lt = 0;
	double Wt = 0;

	Lt = out_p->OmegaT*PI*in_p.R*in_p.R*in_p.D/in_p.B;
	Wt = sqrt(PI*in_p.R*in_p.R); 
	out_p->E_r = Wt/Lt;

}

