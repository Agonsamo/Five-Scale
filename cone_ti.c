/*******************************************************************/
/*  cone_ti.c ->void CONE_TI(in,out)                               */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the sunlit area of one cone that    */
/*  Latest update                              April 20, 1999      */
/*******************************************************************/
/* Output: out_p->tic and out_p->tib                               */
/* A short explanation of this subroutine is found in              */
/* Appendix A of Chen and Leblanc, 1997.	                       */
/*******************************************************************/

# include <math.h>
# include <stdio.h>
# include "data.h" 

void CONE_TI(in_p,out_p)
 
struct PARAMETER in_p;
struct RESULT *out_p;


{  
    double sqrt();
    double EQUATION1();
    double TRIANGLE();

    double m,a,b,c;
    double xa=0,ya=0,xb=0,xd=0,yd=0,xe=0,ye=0,xf=0,yf=0,xg=0,xf2=0,yf2=0; 
    double gamma=PI/2.;
    double B=0,C=0,xe2=0,ye2=0;
    double m1,m2,b1,b2=0;
    double A1=0.,A2=0.,A3=0.,A4=0.; 
    double sza=0,phi=0,vza=0,alpha=0,Hb=0,r=0;
	double arg1=0;
	int cas =0;

	sza=in_p.SZA;
	alpha=in_p.ALPHA;
	vza=out_p->vza;
	phi=out_p->phi;
	Hb=in_p.Hb;
	r=in_p.R;
	

    if ((sza ==0.) && (phi ==0.) )  phi = PI;
  
	
    if(alpha<=sza) 
	{
		arg1 = tan(alpha)/tan(sza);
		gamma = asin(arg1);
	}

    xa=r*cos(vza); ya=0;      /*  2nd radius of the ellipse */

    yf = r*cos(gamma-phi);/* this is the y component of upper part of shadow */

    xf = xa*sin(gamma-phi );

    xb=r*sin(vza)/tan(alpha);

    xg= xb;


    out_p->tib =  2*r*sin(vza)*Hb*(1-phi/PI);




  /***************************** case 1  2  3 *******************************/
	if (vza == 0)         
    {       
		out_p->tic = (PI/2. + gamma)*r*xa ;
		cas=123;
    }
 
 /******************************** 4 5  ************************************/ 
	if (sza >=0 && sza < alpha && vza > 0 && vza < alpha)  
	{       
		out_p->tic =  PI*r*xa ;   
		cas=45;
	}
 /********************************  6  *************************************/
	if (sza >= alpha && sza < PI/2 && vza >0 && vza <= alpha)  
	{       
		cas =6;
        yf2=-r*cos(gamma +phi);
        xf2= xa*sin(gamma+phi);

        yf=r*cos(gamma -phi);
        xf=xa*sin(gamma-phi);

       xg=xb;
                  /*   A* ... are 4 quadrants of ellipse ... */

        if (yf2< 0) 
		{
			A1 = -xa/(2.*r)*EQUATION1(yf2,r);
			A1 += yf2*yf2/(2.*yf2/xf2) ; 
		}

        if (xf>0 ) 
		{
			A2 = xa/(2.*r)*EQUATION1(yf,r);
			A2 -= yf*yf/(2.*yf/xf);   
		}
           
        if(xf>0. && yf2>0.)
		{
			A2 -= xa/(2.*r)*EQUATION1(yf2,r);
			A2 += yf2*yf2/(2.*yf2/xf2) ;
		}

        if (xf <=0. && yf2>0.)
		{
			A2  = PI*r*xa/4.;
			A2 -= xa/(2.*r)*EQUATION1(yf2,r);
			A2 += yf2*yf2/(2.*yf2/xf2)  ; 
		}
        if(xf<=0. && yf2 <=0) A2 = PI*r*xa/4. ;
        if(xf<=0 && xf2<=0) A2 =0;


        if(xf<=0 && yf>0) 
		{
			A3 = -r/(2.*xa)*EQUATION1(xf,xa);
			A3 += yf/xf*xf*xf/2. ; 
		}
   
        if(xf<=0. && yf>0 && xf2 <0.) 
		{
			A3 -= r/(2.*xa)*EQUATION1(xf2,xa);
			A3 += yf2/xf2*xf2*xf2/2.;
		}

        if(xf<=0. && yf<0 && xf2 <0.)
		{
			A3  = PI*r*xa/4.;
			A3 += r/(2.*xa)*EQUATION1(xf2,xa);
			A3 -= yf2/xf2*xf2*xf2/2.;
		}

     
        if(xf<=0. && yf <=0 && xf2 > 0. ) A3 =   PI*r*xa/4. ;




        if (xf <=0 && yf< 0) 
		{
			A4 =  -xa/(2.*r)*EQUATION1(yf,r);
			A4 -=yf*yf/(2.*yf/xf) ;   
		}

       
        if(yf>0) B= TRIANGLE(xf,yf,xg,0.,0.,0.) ;
        if(yf2<0)C= TRIANGLE(xf2,yf2,xg,0.,0.,0.);  
        if(yf2>0)C= - TRIANGLE(xf2,yf2,xg,0.,0.,0.);
        if(yf<0) B= - TRIANGLE(xf,yf,xg,0.,0.,0.) ;
              
		if(yf>-0.00000000001 || yf<0.000000001) B =0; /* added April 20, 1999 to solve a bug when sza=alpha */
		if(yf2>-0.00000000001 || yf2<0.000000001) C =0;
       

        out_p->tic = PI*r*xa;
		out_p->tic -=A1;
		out_p->tic -=A2;
		out_p->tic -=A3;
		out_p->tic -=A4;
		out_p->tic +=B;
		out_p->tic +=C;  


	} 
/***************************   7  8  *************************************/

    if ( sza <= alpha && vza > alpha && vza <= PI/2. )  
    {       
		cas =78;
        yd=r*(1-2*xa*xa/(xb*xb+xa*xa));

        out_p->tic =  PI*r*xa;
		out_p->tic += 2*xb/r*(r*yd-yd*yd/2);
		out_p->tic -= xa/(r)*EQUATION1(yd,r);      
    } 
  
   
 /*******************************   9  **********************************/
    if (sza >= alpha && sza < PI/2 && vza > alpha && vza <= PI/2) 
    { 
		cas =9;
     
        yd=r*(1-2*xa*xa/(xb*xb+xa*xa));
        xd=2*xa*xa*xb/(xb*xb+xa*xa);

        yf2= -r*cos(gamma +phi);
        xf2=  xa*sin(gamma+phi);

        xg=xb;
       
        
        m1 = yf/(xf-xg);
        b1 = -m1*xg;

        a=(r*r+xa*xa*m1*m1);
        b=2.*xa*xa*m1*b1;
        c=xa*xa*(b1*b1-r*r);

        if( b*b>4.*a*c)  xe= (-b + sqrt( b*b-4.*a*c))/(2.*a);
        if( b*b<=4.*a*c)  xe= -b/(2.*a);
        ye= m1*xe+b1;


        m2 = (yf2)/(xf2-xg);
        b2 = yf2-m2*xf2;

        a=(r*r+xa*xa*m2*m2);
        b=2.*xa*xa*m2*b2;
        c=xa*xa*(b2*b2-r*r);
 
        if( b*b>4.*a*c)  xe2= (-b + sqrt( b*b-4.*a*c))/(2.*a);
        if( b*b<=4.*a*c) xe2 = -b/(2.*a);
        ye2= m2*xe2+b2;


        if (yf2 <0) A1 = xa/(2.*r)*(EQUATION1(ye2,r)-EQUATION1(yf2,r)) 
         -(ye2*ye2-yf2*yf2)/(2.*m2) - (-b2*ye2/m2 + b2*yf2/m2);
        if (yf2>=0) A1 =0;

        if(xf > 0  && yf>0) A2 = r/(2.*xa)*EQUATION1(xe,xa) 
         -m1*xe*xe/2. -b1*xe -r/(2.*xa)*EQUATION1(xf,xa) 
         + m1*xf*xf/2. +b1*xf ; 

        if(xf<0 && yf>0)  A2 =  r/(2.*xa)*EQUATION1(xe,xa) 
         -m1*xe*xe/2. -b1*xe  ;

        if(xf<0 && yf <=0) A2 = PI*r*xa/4.; 


        if(yf2>0 && xf2>0) A2 = A2 -  xa/(2.*r)*(EQUATION1(yf2,r)-EQUATION1(ye2,r)) 
         +(yf2*yf2-ye2*ye2)/(2.*m2) + (-b2*yf2/m2 + b2*ye2/m2);
         
        if(yf2>0 && xf2<=0)  A2 = A2 -  r/(2.*xa)*(EQUATION1(xe2,xa)) 
         +( m2*xe2*xe2/2. + b2*xe2) ; 

 
        if(xf<0 && yf>0)  A3 =  -r/(2.*xa)*EQUATION1(xf,xa) +m1*xf*xf/2. +b1*xf;

        if(xf<0 && yf <=0) A3 = PI*r*xa/4.;

        if(yf2>0 && xf2<=0) A3 = A3 - ( - r/(2.*xa)*EQUATION1(xf2,xa)
                                      +m2*xf2*xf2/2. +b2*xf2) ;


    
        if(yf<0 && xf<0) A4 =    -xa/(2.*r)*EQUATION1(yf,r)
         -(yf*yf)/(2.*m1) + (b1*yf/m1) +(b1*b1)/(2.*m1) - (b1*b1/m1) ;   

        if(yf<0 && xf <0) A1 = -xa/(2.*r)*EQUATION1(ye,r) + xe*(ye-b1)/2. ; 
              
       
         
        C = 0 ;
        if (xf< 0  && xf2 > 0 ) 
          {

           m= (yd-ye)/(xd-xe);
           b= yd -m*xd;
               C = TRIANGLE(xe,ye,xd,yd,xg,0.) - 
               xa/(2.*r)*(EQUATION1(yd,r)-EQUATION1(ye,r))   
               +(yd*yd/2.-yd*b)/m - (ye*ye/2.-ye*b)/m ;


           }

			if ( xf <0 && xf2 <= 0 && yf2>0)

			{
			if ((xe2-xe)< 0.0000001 &&  (xe2-xe) > -0.00000001)
			{
				m=0;
				b= ye; 
				C= TRIANGLE(xe,ye,xe2,ye2,xg,0.)-
				2*(xa/(2.*r)*EQUATION1(ye2,r) - xe*(ye2)) ;
			}else
			{
				m= (ye2-ye)/(xe2-xe);
				b= ye -m*xe;
				C = TRIANGLE(xe,ye,xe2,ye2,xg,0.) - 
					xa/(2.*r)*(EQUATION1(ye2,r)-EQUATION1(ye,r))   
					+(ye2*ye2/2.-ye2*b)/m - (ye*ye/2.-ye*b)/m ;

			} 

           
		}

		out_p->tic = PI*r*xa;
		out_p->tic +=2*xb*(r*yd-yd*yd/2.)/r;
		out_p->tic -= xa/r*EQUATION1(yd,r);
		out_p->tic -=A1;
		out_p->tic -=A2;
		out_p->tic -=A3;
		out_p->tic -=A4;
		out_p->tic -=C ; /* because of some strage problems with VC++, I had to put the tic computation  on 8 different lines */ 
	


	}                            /* end of  case 9 */

   
 /*************************************************************************/


}                               /* end of all !! */
     

