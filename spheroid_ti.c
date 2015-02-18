/*******************************************************************/
/*  spheroid_ti.c ->void SPHEROID_TI(in,out)                       */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates sunlit crown proportion             */
/*  for deciduous (spheroid) shape (based on Li & Strahler )       */ 
/*  Latest update                                    June 4, 1996  */
/*******************************************************************/


# include <stdio.h>
# include <math.h>
# include "data.h"


void SPHEROID_TI(in_p,out_p)
 
struct PARAMETER in_p;
struct RESULT *out_p;



{  

 double  vza_prime,sza_prime,cs_prime;


         vza_prime = atan(in_p.Hb/(2*in_p.R)*tan(out_p->vza));
         sza_prime = atan(in_p.Hb/(2*in_p.R)*tan(in_p.SZA));

          cs_prime = cos(sza_prime)*cos(vza_prime) + 
                     sin(sza_prime)*sin(vza_prime)*cos(out_p->phi);

           out_p->tib = out_p->tab*0.5*(1.+ cs_prime); 

         


 }                    

      

    
     

