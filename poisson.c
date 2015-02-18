/*******************************************************************/
/*  poisson.c ->void POISSON(in,out)                               */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the poisson distribution for random */
/*  distribution of trees                                          */
/*  Latest update                                February 4, 1997  */
/*******************************************************************/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include "data.h" 

void  POISSON(in_p,out_p)

struct PARAMETER in_p;
struct RESULT *out_p; 


{
	double exp();
	double M;                  /* Mean value of the distribution */ 
	int ai,ak;    		 	                  /* dummy variables */
	double PNN=0;          /* initial value in computaion of *PX */
	double PX_tot=0;                           /* sum of all *PX */ 
	
	M = in_p.D/in_p.n; 

	for (ak=0;ak<NN;ak++)
	{
                             
	   PNN =  exp(-M);
       	   for (ai=1;ai<=ak;ai++)
	   {
			PNN =PNN*M/(ai*1.);
	   }
			out_p->Px[ak] = PNN; 
			PX_tot = PX_tot + out_p->Px[ak]; 
                                 
	}
	if (PX_tot <0.95)
	{ 
		printf("\n Possible error in Poisson distribution:");
		printf(" sum of all tree distribution probabilities is less ");
		printf("than 0.95 (sum=%8.5f) \n",PX_tot);
	} 

	for (ak=0;ak<NN;ak++)  out_p->Px[ak] = out_p->Px[ak]/ PX_tot;


   /* in case IMAX is not large enough so that PX_tot equals unity,
       a normalisation is done for PX */ 

} 
