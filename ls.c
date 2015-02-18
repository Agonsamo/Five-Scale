/*******************************************************************/
/*  ls.c ->void LS(in,out,shape)                                   */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  This subroutine computes Cv, Cs, Lo_90, mu                     */
/*  Latest update                                   January 2004   */
/*******************************************************************/
 
 
# include <stdio.h>
# include <math.h>
# include "data.h"
 
void LS(in_p,out_p,shape)
 
struct PARAMETER in_p;
struct RESULT *out_p;
int shape; 

{
	double LAI_HOT_SPOT =1;  

	out_p->mu = in_p.LAI*in_p.B/(in_p.D*out_p->V);

	out_p->Ls=out_p->Gs*LAI_HOT_SPOT/cos(in_p.SZA);   

	out_p->H = 1./(out_p->mu);   
	

	out_p->lambda_m = out_p->H*tan(out_p->xi) ; 

	if(shape==1)
	{
		out_p->Lo_90 = out_p->mu*PI*in_p.R*(in_p.Hb+out_p->Hc/3.)/(2*in_p.Hb+out_p->Hc);
	

		out_p->Cs = out_p->Gs*out_p->Ss*out_p->mu/out_p->Lo_90;
		out_p->Cv = out_p->Gv*out_p->Sv*out_p->mu/out_p->Lo_90;  /*test 5 Feb 1999 */

	}

	if(shape==2) 
	{
		out_p->Lo_90 = out_p->mu*in_p.R*4/3.; /*  January 12, 2004 */


		out_p->Cs = out_p->Gs*out_p->Ss*out_p->mu/out_p->Lo_90;
		out_p->Cv = out_p->Gv*out_p->Sv*out_p->mu/out_p->Lo_90; 
		
	}

}
