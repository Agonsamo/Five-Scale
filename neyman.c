/*******************************************************************/
/*  neyman.c ->void NEYMAN(in,out)                                 */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  This subroutine computes the Neyman distribution used to       */
/*  simulate clumping  of trees                                    */
/*  Latest update                              February 12, 1998   */
/*******************************************************************/



# include <math.h>
# include <stdio.h>
# include "data.h" 


void NEYMAN(in_p,out_p)  

struct PARAMETER in_p;
struct RESULT *out_p; 

{  

  double pow(),exp();
  double m1,Px_tot=0;
  double PXX;
  int t;

  register int k,i;

 

  m1 = in_p.D/(in_p.n*in_p.m2); 

  out_p->Px[0] = exp(-m1*(1.-exp(-in_p.m2))); 
  
  Px_tot=  out_p->Px[0]   ;

  for (k=1;k<NN;k++)
  {
    out_p->Px[k] = 0;

 /*********************************************************************/
 /*  to increse the speed, the model                                  */
 /*  does not compute values smaller than 0.00001 for probabilities   */
 /*  greater than the mean number of tree stand per quadrat           */
 /*********************************************************************/

    if(k> in_p.D/in_p.n && out_p->Px[k-1] < 0.00001) 
    {
		out_p->Px[k] =0; 
    } else
    {
		for(t=0;t<k;t++)
        { 
			PXX =    m1*in_p.m2*exp(-in_p.m2)/(k)*out_p->Px[k-t-1] ;               
			for(i=1;i<=t;i++)
            {
				PXX = PXX*in_p.m2/i*1.;                  
            }

			out_p->Px[k] = out_p->Px[k] + PXX;

        }  
	}

 /*********************************************************************/
 /*************** making sure that the sum of Px is unity ***************/
 /*********************************************************************/
    Px_tot = Px_tot + out_p->Px[k];
 
  }

if(Px_tot < 0.98)
{
 printf("\n Possible error in Neyman distribution Px_tot = %8.5f\n",Px_tot);
 out_p->Error[1] =1;
 
}
  for (k=0;k<NN;k++)
   {
       out_p->Px[k] = out_p->Px[k]/Px_tot;
   } 

} 
