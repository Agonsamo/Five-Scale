/*******************************************************************/
/*  cone_ta.c ->void CONE_TA(in,out)                               */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates crown projection to viewer and the  */
/*  crown volume                                                   */
/*  Latest update                              February  5, 1997   */
/*******************************************************************/
/*   3 possible cases:                                             */
/*   1 --> nadir (VZA =0) 					                       */
/*   2 --> half apex angle > VZA > 0                               */
/*   3 --> VZA >= half apex angle                                  */
/*******************************************************************/


# include <math.h>
# include "data.h" 
# include <stdio.h>

void CONE_TA(in_p,out_p)

struct PARAMETER in_p;
struct RESULT  *out_p;


{  

	double EQUATION1();
	double gamma=0.,xa=0.,xb=0.,yd=0. ;
   
  

	if(in_p.ALPHA <= out_p->vza) gamma= asin(tan(in_p.ALPHA)/tan(out_p->vza));
   
 /****************************** 1  ************************************/
    if (out_p->vza ==0)         
    {   
      out_p->tac = PI*in_p.R*in_p.R;   
      out_p->tab = 0 ;  
    }
 /******************************** 2  ***********************************/

    else if (out_p->vza > 0 && out_p->vza <= in_p.ALPHA)   
    {   
      out_p->tac = PI*in_p.R*in_p.R*cos(out_p->vza) ;
      out_p->tab = 2*in_p.R*in_p.Hb*sin(out_p->vza); 
	
    }

	
 /******************************** 3  ***********************************/

	else if (out_p->vza >  in_p.ALPHA ) 
    {
		xa=in_p.R*cos(out_p->vza);  
		xb=in_p.R*sin(out_p->vza)/tan(in_p.ALPHA);
		yd=in_p.R*(1-2*xa*xa/(xb*xb+xa*xa));

		out_p->tac = PI*in_p.R*xa;
		out_p->tac += 2*xb/in_p.R*(in_p.R*yd-yd*yd/2);
		out_p->tac -= xa/(in_p.R)*EQUATION1(yd,in_p.R);

		out_p->tab = 2*sin(out_p->vza)*in_p.R*in_p.Hb;
             
       
    }     
    
 /************************** crown volume *******************************/

 out_p->V= PI*in_p.R*in_p.R*(in_p.Hb+out_p->Hc/3.); /* January 12, 2004 */
}    


    
