/*******************************************************************/
/*  equation.c ->double EQUATION1(x,xa)                            */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Used in the solution for an integral                           */
/*  Latest update                                September 1998    */
/*******************************************************************/



# include <math.h>
# include <stdio.h>

double  EQUATION1(x,xa)
double x,xa;
 
{
  
   double sqrt(),asin();
   double answer=0;  


   if( xa*xa>x*x) answer = x*sqrt(xa*xa-x*x) + xa*xa*asin(x/xa);
   if (xa*xa<=x*x) answer =0; 



   return answer; 
  
 }
