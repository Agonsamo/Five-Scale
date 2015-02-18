/*******************************************************************/
/*  triangle.c ->double  TRIANGLE(xa,ya,xb,yb,xc,yc)               */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the area of any triangle            */
/*  Latest update                              November 27, 1995   */
/*******************************************************************/



# include <math.h>
# include <stdio.h>


 double  TRIANGLE(xa,ya,xb,yb,xc,yc) 
 double xa,xb,xc,ya,yb,yc;  /* coordinates of the 3 corners of the triangle */
 
{
  
   double pow(),sqrt(),acos();
   double a,b,c,C,cC,answer;  

   a = sqrt((yc-yb)*(yc-yb)+(xc-xb)*(xc-xb));
   b = sqrt((yc-ya)*(yc-ya)+(xc-xa)*(xc-xa));
   c = sqrt((yb-ya)*(yb-ya)+(xb-xa)*(xb-xa));

     cC= (a*a+b*b-c*c)/(2.*a*b); 
      C= acos(cC);


   answer = 0.5*a*b*sin(C);

   return answer; 
  
 }
