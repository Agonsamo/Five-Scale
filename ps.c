/*******************************************************************/
/*  ps.c ->void PS(in,out)                                         */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates some parameters for hotspot kernel  */ 
/*  Latest update                                      June, 1997  */
/*******************************************************************/

 
 
# include <stdio.h>
# include <math.h>
# include "data.h"
 
void PS(in_p,out_p)

struct PARAMETER in_p;
struct RESULT *out_p;
 
{

	out_p->PS = out_p->Pig*out_p->Pvg;
	if (out_p->Sg_0 <= 0) out_p->Wt = 0;
	else  out_p->Wt = sqrt(out_p->Sg_0);

	if (in_p.SZA == PI) out_p->H=0 ;  
	if (in_p.SZA != PI) out_p->H=1./cos(in_p.SZA)*(out_p->Hc/3. + in_p.Hb + in_p.Ha);
          
 
	out_p->Lt= 1.*in_p.D*out_p->Sg_0/in_p.B*out_p->OmegaT;   
	if (out_p->xi < PI)  out_p->lambda_m = out_p->H*tan(out_p->xi);      
	if (out_p->xi >=PI)  out_p->lambda_m =0;
       
    
}

      

