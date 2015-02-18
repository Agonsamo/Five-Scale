/*******************************************************************/
/*  q1_ms.c ->void Q1(in,out,DEL)                                  */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that calculates appr. of Q1tot for MS scheme        */    
/*  Latest update                              September 23, 1998  */
/*******************************************************************/

 
# include <math.h>
# include <stdio.h>
# include "data.h"


double Q1(theta,out_p,in_p,phi)
 
struct PARAMETER in_p;
struct RESULT *out_p;
double theta;
double phi;


{

	double Cv=0;
	double Q1_theta=0;
	double PgapV_theta=0;
	double Vg_theta=0;
	double DEL = 0;
	double diff=0;


	if(theta<0)theta=-theta;
	diff = cos(in_p.SZA)*cos(theta) +sin(theta)*cos(phi)*sin(in_p.SZA);/*XI(in_p.SZA,theta,phi);*/
	if (diff<0) diff = - diff;

	DEL = 1-diff*in_p.Cp/PI;  
	if(DEL<0) DEL=-DEL;
	if(DEL>1) DEL=1;

	if(theta>PI/2) theta=PI-theta; 

	
	if(theta==0) theta = 0.0000000001; /* 25 July 1999, solves some of the warning 14 problems */
	if(sqrt((theta-PI/2)*(theta-PI/2))<0.000000001) theta=PI/2-0.0000001; 
	Cv =out_p->Gv/sin(theta);
 
	Vg_theta = 2*tan(theta)*in_p.R*(in_p.Hb+out_p->Hc) +PI*in_p.R*in_p.R;

	PgapV_theta = exp(-Cv*in_p.LAI*in_p.B/(in_p.D*Vg_theta*cos(theta)));

	

	Q1_theta=(1.-exp(-(out_p->Cs*out_p->Lo_90+Cv*out_p->Lo_90))) *out_p->Cs*Cv/(out_p->Cs+Cv)/(1-PgapV_theta);
	
	Q1_theta = Q1_theta*DEL; 

/*	fprintf(fp,"\nVZA: %6.1f d=%7.2f THETA= %7.3f Q1_theta: %6.3f DEL: %6.3f  rtheta %6.3f diff %6.3f phi %6.3f",
		out_p->vza*180/PI,out_p->E_r,theta*180/PI,Q1_theta,DEL*180/PI,real_theta*180/PI,diff*180/PI,phi*180/PI);
	fclose(fp); */

	out_p->Error[14] =0;
	out_p->Error[16] =0;

	if(Q1_theta >1) 
	{
		out_p->Error[16] =1;
		return 1;
	} 
	else if(Q1_theta <0) 
	{	
		out_p->Error[14] =1;

		return 0;	
	} 
	else return Q1_theta;


}
 
 