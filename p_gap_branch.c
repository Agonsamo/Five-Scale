/*******************************************************************/
/*  p_gap_branch.c ->void P_GAP_BRANCH(in,out,choice)              */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates the gap fraction in one crown       */
/*  Latest update                            December 2, 1998      */
/*******************************************************************/ 
/*    ZA     = solar or view zenith angle 						   */
/*    in_p.B      = size of the domain					           */
/*    in_p.D      = number of tree stand in the domain			   */
/*    V      = volume of one crown					               */
/*    s      = s(VA) mean path through a crown				       */ 
/*    in_p.Ll     = branch leaf area index    	                   */
/*    Lb     = branch silhouette area index				           */
/*    muB, muL = density of Lb and in_p.Ll inside crown            */
/*    in_p.ALPHA_B = branch angle to the horizontal			       */
/*    in_p.ALPHA_L = leaf angle to the horizontal			       */
/*    in_p.R = radius of the crown (minor axis if spheroid shape)  */
/*    in_p.GAMMA_E = needle-to-shoot   				               */
/*    in_p.RATIO  = leaf thickness to width-length RATIO		   */ 
/*******************************************************************/


# include <stdio.h>
# include <math.h>
# include <string.h>
# include "data.h"

void P_GAP_BRANCH(in_p,out_p,CHOICE)

struct PARAMETER in_p;
struct RESULT *out_p;
char CHOICE[];
{
	double exp();
 
	double xb=0,xl=0,Lb=0,ZA=0;   
	double Gb=0,Gl=0,Pbj=0,Pl1=0,Pj=0;  
	double cos_ZAb=0,P=0,mub=0,mul=0,beta=0,s=0;
	int j=0;
	double aGb=0, aGl=0;

	Lb =in_p.LAI/in_p.Ll;                   /* Eq. 8 [2] */ 
	mub=Lb*in_p.B/(out_p->V*in_p.D);        /* Eq. 7 [2] */  
	mul=in_p.Ll*in_p.B/(out_p->V*in_p.D);

	if(!strcmp("VZA",CHOICE))
	{
		ZA = out_p->vza;
		s = out_p->Sv;
	}
	if(!strcmp("SZA",CHOICE))
	{
		ZA = in_p.SZA; 
		s = out_p->Ss;
	}
	if(!strcmp("0",CHOICE))
	{
		ZA = 0; 
		s = 1/3*out_p->Hc+in_p.Hb;
	}

	if(!strcmp("LAI",CHOICE))
	{
		ZA =acos(0.537+0.025*in_p.LAI);  
		s = out_p->V/(out_p->Vg_0_mean*(0.537+0.025*in_p.LAI));
	}





	if(in_p.ALPHA_L>0)
	{
		if (ZA<=PI/2.-in_p.ALPHA_L) Gl = cos(in_p.ALPHA_L)*cos(ZA) ;

		if (ZA>PI/2. -in_p.ALPHA_L)
		{
			xl = 1/tan(in_p.ALPHA_L)*1/tan(ZA);
			xl = acos(xl);
			Gl = cos(in_p.ALPHA_L)*cos(ZA)*(1.+2.*(tan(xl)-xl)/PI) ;
		}

		aGl = acos(Gl); 
		Gl=Gl + sin(aGl)*in_p.RATIO;

	} else
	{

		Gl = 0.5;
		aGl =acos(Gl);
	}

	if(!strcmp("VZA",CHOICE))
	{
		out_p->GFoliage = Gl;
	}

	if(in_p.ALPHA_B>0)
	{
             /***** Gb computation  Eq. 9-11 [2]   *****/   

		if (ZA<=PI/2.-in_p.ALPHA_B) Gb = cos(in_p.ALPHA_B)*cos(ZA) ;

		if (ZA>PI/2. -in_p.ALPHA_B)
		{
			xb = 1/tan(in_p.ALPHA_B)*1/tan(ZA);
			xb = acos(xb);
			Gb = cos(in_p.ALPHA_B)*cos(ZA)*(1.+2.*(tan(xb)-xb)/PI) ;
		}

		aGb = acos(Gb);
		Gb=Gb + sin(aGb)*cos(in_p.ALPHA_B)*in_p.Rb/in_p.R; /* Eq. 11  */

	}

	else Gb =0.5;

    Pl1 = 0;
    for (beta=0;beta<=PI;beta=beta + PI/100.)
    {
		/*  Equation 12  [2]    */ 
		if(in_p.ALPHA_B>0)
		{
			cos_ZAb = sin(ZA)*sin(in_p.ALPHA_B)*cos(beta)+cos(ZA)*cos(in_p.ALPHA_B);
			aGb = acos(cos_ZAb);
			cos_ZAb = cos_ZAb + sin(aGb)*cos(in_p.ALPHA_B)*in_p.Rb/in_p.R;
		}
		else cos_ZAb =1;

		if (cos_ZAb<0) cos_ZAb = -cos_ZAb;

		/*  Equation 13  [2]    */ 
		aGl = Gl*in_p.Ll/(in_p.GAMMA_E*cos_ZAb);
		Pl1 = Pl1 + 1/PI*exp(-aGl)*PI/100.; 

    }
    
    Pbj=exp(-Gb*mub*s*1.); 
    P=Pbj; 

    for (j=1;j<=s*5.;j++)   /* This limits the number of branches that */
                            /* can be encountered  */
    {

		Pbj= Pbj*Gb*mub*s*1./j;
		Pj = Pbj*pow(Pl1,j*1.);

		if (Pj <= 1.) P = P + Pj;   /**** Equ. 17  [2]   ****/  
    }

	if(!strcmp("VZA",CHOICE))
    {
		out_p->PgapV=P;
		out_p->Gv= cos(ZA)*out_p->Vg_0*in_p.D*log(1/P)/(in_p.B*in_p.LAI); 
    }

	if(!strcmp("SZA",CHOICE)) 
    {
		out_p->PgapS=P;
		out_p->Gs= cos(ZA)*out_p->Sg_0*in_p.D*log(1/P)/(in_p.B*in_p.LAI); 
		out_p->Lo = in_p.LAI*in_p.B/(in_p.D*out_p->Sg_0*cos(ZA));
    }
	if(!strcmp("0",CHOICE)) 
    {
		out_p->Pgap0=P;
    }

	 if(!strcmp("LAI",CHOICE)) 
    {
		out_p->PgapV_mean=P;
    }

}  
