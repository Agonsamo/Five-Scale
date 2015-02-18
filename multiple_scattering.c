/*******************************************************************/
/* multiple_scattering.c ->void MULTIPLE_SCATTERING(in,out) */
/* Jing M. Chen, jing.chen@ccrs.nrcan.gc.ca */
/* Sylvain G. Leblanc sylvain.leblanc@ccrs.nrcan.gc.ca */
/*******************************************************************/
/* This subroutine computes the amount of electromagnetic */
/* radiation reaching the four components (sunlit and shaded */
/* background and foliage) due to multiple scattering */
/* It serve as the basis for the hyperspectral mode of 5-Scale */
/* Latest update March 2001 */
/*******************************************************************/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>
# include "data.h"
void MULTIPLE_SCATTERING(in_p,out_p)
struct PARAMETER in_p;
struct RESULT *out_p;
{
double pow(),exp(),atan(),cos(),sin(),tan();
double Q1(); /* amount (%) of sunlit foliage seen in one crown */
double XI();
double DEL=0;
double DEL2=0;
int i=0;

int wave_index_min =1; /* default band used , Red, NIR, SWIR) */
int wave_index_max =4;
double TG =0; /* to check tan sign */
double Intensity=0;
double argt=0;
double arg3=0, arg4=0, arg5=0; /* used in different exp() */
double dem=0;
double CV=0; /* similar to Cv in [1] (old definition) */
double Cv_mean=0;
double F_s_TRest =0; /* sky view factor to
foliage not weighted with height */
double F_s_T=0; /* sky view factor to ALL foliage */
double F_T_T=0; /* total foliage view factor from foliage */
double F_G_T=0; /* total background view factor from foliage */
double F_g_T=0; /* sunlit background view factor from foliage */
double F_zg_T=0; /* shaded background view factor from foliage */
double F_t_t=0; /* sunlit foliage view factor from sunlit foliage */
double F_zt_t=0; /* used for other to find F_t_t */
double F_t_zt=0; /* sunlit foliage view factor from shaded foliage */
double F_s_G=0; /* sky view factor from ground */
double F_t_G=0; /*sunlit foliage view factor from ground */
double F_tt2=0;
double F_tt3=0;
double F_tt_zt=0;/*shaded side of sunlit leaves v.f. from shaded foliage */

double F_tt_G=0; /*shaded side of sunlit leaves v.f. from background */
int p; /*index for Phi ...*/
double fd[MAX_WAVE];
double GZ[MAX_WAVE];
double Hi=0,Hj=0;
double lambda[MAX_WAVE];
double Li=0, Lj =0; /* Leaf area index based on height Hi, Hj ...*/
double num=0;
double OMEGA_TOTAL =0;
double PgapV_mean=0;
double P_tot=0;
double Q1_mean=0,Q2_mean=0;
double Q1B_mean=0,Q2B_mean=0;
double Rr=0;
double RG_2nd[MAX_WAVE]; /* second order of scattered light reflected by
the ground */
double RT_2nd[MAX_WAVE];
double RTZ_2nd[MAX_WAVE];
double RGT_2nd[MAX_WAVE];
double RGT_3rd[MAX_WAVE];
double RT_MS[MAX_WAVE];
double RG_MS[MAX_WAVE];
double RZT_MS[MAX_WAVE];
double RT_HC[MAX_WAVE];
double RG_HC[MAX_WAVE];
double RZT_HC[MAX_WAVE];
double GZ_HC[MAX_WAVE];
double TT_HC[MAX_WAVE];
double M_factor[MAX_WAVE];
/*Added by Ting, a new M_factor*/
double M_factor2[MAX_WAVE];
/**/
double TT[MAX_WAVE];
double TZ[MAX_WAVE];
double theta=0;
double theta_min=0;
double theta_max=0;
double theta_half_mean = 0;
double theta_h;
double inc_theta;
double tauR =0;
double Vg_mean=0;


/***** init of shaded reflectivities in the out_p structure **/

out_p->TZ1 =0 ;
out_p->TZ2 =0 ;
out_p->TZ3 =0 ;
out_p->TZ4 =0 ;
out_p->GZ1 =0 ;
out_p->GZ2 =0 ;
out_p->GZ3 =0 ;
out_p->GZ4 =0 ;

Intensity = atan((in_p.Hb + out_p->Hc)/(2*in_p.R));
Intensity = Intensity-in_p.SZA;
if(Intensity<0) Intensity=-Intensity;
Intensity = cos(Intensity);

/********************** Sky diffusion **************************/

if(in_p.SPECTRAL>1)
{
wave_index_min =0;
wave_index_max =152;



} else
{
wave_index_min =1;
wave_index_max =5;
in_p.RG[1] = in_p.G1;
in_p.RG[2] = in_p.G2;
in_p.RG[3] = in_p.G3;
in_p.RG[4] = in_p.G4;

in_p.RT[1] = in_p.T1;
in_p.RT[2] = in_p.T2;
in_p.RT[3] = in_p.T3;
in_p.RT[4] = in_p.T4;

in_p.TT[1] = in_p.TT1;
in_p.TT[2] = in_p.TT2;
in_p.TT[3] = in_p.TT3;
in_p.TT[4] = in_p.TT4;

}



for(i=wave_index_min;i<wave_index_max;i++)
{
if(in_p.SPECTRAL>1)
{
lambda[i] = out_p->Wave[i]/1000.;
in_p.RG[i] = out_p->BACKGROUND_REF[i];
in_p.RT[i] = out_p->FOLIAGE_REF[i];
in_p.TT[i] = out_p->FOLIAGE_TRANS[i];
}
else
{
lambda[i] = in_p.Band[i]/1000.;
}


/* rayleigh Scattering */
tauR = 0.008569*pow(lambda[i],-4)*(1-0.0113*pow(lambda[i],-2.)+
0.00013*pow(lambda[i],-4));
Rr = (3*tauR+(2-3*cos(in_p.SZA))*(1-exp(-tauR/cos(in_p.SZA))))/(4+3*tauR);
fd[i] = ((1-Rr) - exp(-tauR/cos(in_p.SZA)))/(1-Rr);

//Added By Ting to adjust the aerosol effects on the diffuse fraction: 2 times or 3times.
fd[i] = fd[i]*2;
//End
// Added by Ting to simulate the MODIS senario
//fd[i] = 0;
//End

}

if(in_p.LAI>0)
{

/********************** OMEGA TOTAL ***********************/
num = log(out_p->Pvg_mean);
dem = log(exp(-0.5*in_p.LAI/(0.537+0.025*in_p.LAI)));
OMEGA_TOTAL = num/dem;

out_p->Error[6] =0; /* Win 95/NT Warnings */
out_p->Error[5] =0;
if(OMEGA_TOTAL > 1)
{
OMEGA_TOTAL = 1;
out_p->Error[6] =1;
}
CV = 0.5*OMEGA_TOTAL/cos(out_p->vza);
/**************************** Q1 and Q2
**************************************/
theta_half_mean =
atan((in_p.Ha+0.5*(in_p.Hb+out_p->Hc/3))/(0.5*out_p->E_r));
Q1_mean = Q1(PI/2+theta_half_mean,out_p,in_p,0);


DEL = 1-in_p.Cp*XI(in_p.SZA,PI/2+theta_half_mean,0)/PI;
Q1B_mean = Q1_mean*(1-DEL)/DEL;
if(Q1B_mean >1) Q1B_mean =1;
if(Q1B_mean <0) Q1B_mean =0;

Cv_mean = out_p->Gv*out_p->Sv*out_p->mu/out_p->Lo_90;
TG = tan(PI/2+theta_half_mean);
if (TG < 0) TG = -TG;
Vg_mean = 2*TG*in_p.R*(in_p.Hb+out_p->Hc/3) + PI*in_p.R*in_p.R;
if(Vg_mean < 0) Vg_mean =0;
PgapV_mean =
exp(-Cv_mean*in_p.LAI*in_p.B/(in_p.D*Vg_mean*cos(theta_half_mean)));

DEL = 1-in_p.Cp*(cos(in_p.SZA)*cos(theta_half_mean)
+sin(theta_half_mean)*sin(in_p.SZA))/PI;
if(DEL<0) DEL=-DEL;
Q2_mean =
out_p->Cs*Cv_mean/(Cv_mean-out_p->Cs)*(exp(-out_p->Cs*out_p->Lo_90)-exp(-Cv_mean*out_p->Lo_90))/(1-PgapV_mean);

Q2B_mean = (1-DEL)*Q2_mean;
Q2_mean = Q2_mean*DEL;
if(Q2B_mean >1) Q2B_mean =1;
if(Q2B_mean <0) Q2B_mean =0;
/*****************************transmitance ***************************/

for(i=wave_index_min;i<wave_index_max;i++)
{
if(strcmp(in_p.SHAPE,"SPHEROID")) TT[i] = pow(in_p.TT[i],in_p.GAMMA_E) ;
else TT[i] = in_p.TT[i];

}





/***************************************************************************************/

F_G_T =0;
for(Li=in_p.DeltaLAI;Li<=in_p.LAI;Li=Li+in_p.DeltaLAI)
{
/***** sky view factor from foliage, weighted by view penetration ******/
arg4 = 0.5*OMEGA_TOTAL*Li/cos(asin(0.537+0.025*Li));
F_s_T = F_s_T + 0.5*exp(-arg4)*exp(-CV*Li)*in_p.DeltaLAI ;
F_s_TRest = F_s_TRest + 0.5*exp(-arg4)*in_p.DeltaLAI;
arg5 =0.5*OMEGA_TOTAL*(in_p.LAI-Li)/cos(asin(0.537+0.025*(in_p.LAI-Li)));
F_G_T = F_G_T + 0.5*exp(-arg5)*in_p.DeltaLAI;
P_tot = P_tot + exp(-CV*Li)*in_p.DeltaLAI;
F_tt2 =0;
F_tt3 =0;
Hi = (in_p.Hb+out_p->Hc/3)-Li*(in_p.Hb+out_p->Hc/3)/in_p.LAI;
inc_theta =in_p.DeltaLAI;
theta_min = (Hi-inc_theta)/out_p->E_r;

theta_min = 0.5*(PI/2+atan(theta_min));

theta_max = ((in_p.Hb+out_p->Hc/3)-Hi)/out_p->E_r;

theta_max = 0.5*(PI/2+atan(theta_max));



if(theta_min>theta_max)
{
theta = theta_min;
theta_min = theta_max;
theta_max =theta;
}
for (theta=theta_min; theta<=theta_max;theta=theta+inc_theta)
{

theta_h = 2*theta-PI/2;

theta_h=PI/2-theta_h;
F_tt2 =
F_tt2 +
Q1(theta_h,out_p,in_p,0)*cos(theta)*sin(theta)*inc_theta*10*PI/180;
DEL =0;
DEL2=0;
for(p=0;p<180;p=p+10)
{
DEL = DEL + 1-in_p.Cp*(cos(in_p.SZA)*cos(theta_h)
+sin(theta_h)*sin(in_p.SZA)*cos(p*PI/180))/PI;
DEL2 = DEL2 + in_p.Cp*(cos(in_p.SZA)*cos(theta_h)
+sin(theta_h)*sin(in_p.SZA)*cos(p*PI/180))/PI;
}
F_tt2 = F_tt2*DEL*10*PI/180;
F_tt3 = F_tt2*DEL2*10*PI/180;

}
F_t_zt = F_t_zt + F_tt2*in_p.DeltaLAI ;
F_tt_zt = F_tt_zt+ F_tt3*in_p.DeltaLAI ;

}


F_t_zt = F_t_zt/(in_p.LAI); /* sunlit foliage view factor from shaded foliage */
F_tt_zt= F_tt_zt/(in_p.LAI); /*shaded side of sunlit leaves (tt means transmit) v.f. from shaded foliage */
F_zt_t = F_t_zt;
F_s_T = F_s_T/P_tot;   //sky view factor for the foliage
F_s_TRest = F_s_TRest/(in_p.LAI); // sky view factor
F_G_T = F_G_T/(in_p.LAI);
F_g_T = F_G_T*out_p->Pig;
F_zg_T = F_G_T*(1-out_p->Pig);
/********* sky view factor from ground ******/
arg3 = 0.5*OMEGA_TOTAL*in_p.LAI/(0.537+0.025*in_p.LAI);
F_s_G = exp(-arg3);   /* sky view factor from ground */
/***** sunlit trees view factor from ground ******/
F_t_G = 0.5*(Q1_mean+Q2_mean)*(1-F_s_G);
F_tt_G = 0.5*(Q1B_mean+Q2B_mean)*(1-F_s_G);

F_t_t = 1-F_s_TRest - F_G_T - F_zt_t;
F_T_T = 1-F_G_T-F_s_TRest;
/* Warning testing (for Windows 95/NT version) */
if(F_t_t < 0)
{
F_t_t =0;
out_p->Error[5] =1;
}
if(F_s_T < 0)
{
F_s_T =0;
out_p->Error[5] =1;
}
if(F_zt_t < 0)
{
F_zt_t =0;
out_p->Error[5] =1;
}
if(F_t_G < 0)
{
F_t_G =0;
out_p->Error[5] =1;
}
if(F_G_T < 0)
{
F_G_T =0;
out_p->Error[5] =1;
}

/******************** calculation of the "shaded reflectivities
***************/
for (i=wave_index_min;i<wave_index_max;i++)
{

RT_2nd[i] = in_p.RT[i]*((in_p.RT[i]+ TT[i]*out_p->Ft)*F_t_t +
in_p.RG[i]*F_g_T + fd[i]*F_s_T);

RTZ_2nd[i] = in_p.RT[i]*(in_p.RT[i]*F_t_zt + TT[i]*F_tt_zt +
in_p.RG[i]*F_g_T + fd[i]*F_s_T);
RG_2nd[i] = in_p.RG[i]*(in_p.RT[i]*F_t_G +TT[i]*F_tt_G + fd[i]*F_s_G);
RGT_2nd[i] = (RG_2nd[i]*(F_G_T) + RTZ_2nd[i]*(1-F_G_T-F_t_t-F_s_TRest) +
RT_2nd[i]*(F_t_t+F_zt_t))/(1.-F_s_TRest);
RGT_3rd[i] =(in_p.RT[i]+TT[i])/2*(1-F_s_TRest)*RGT_2nd[i];
RT_MS[i] = RT_2nd[i] + RGT_3rd[i]/(1-(in_p.RT[i]+TT[i])/2*(1-F_s_TRest));
RZT_MS[i] = RTZ_2nd[i] +RGT_3rd[i]/(1-(in_p.RT[i]+TT[i])/2*(1-F_s_TRest));
RG_MS[i] = RG_2nd[i] +
in_p.RG[i]*RGT_2nd[i]*(1-F_s_G)/(1-(in_p.RT[i]+TT[i])/2*(1-F_s_TRest));

TZ[i] = RZT_MS[i] ;
GZ[i] = RG_MS[i]; /*+ TCrown[i]*in_p.RG[i]; */
/*
out_p->RT_HC[i] = (in_p.RT[i]*Intensity + RT_MS[i]);
out_p->RG_HC[i] = (in_p.RG[i]*(1-out_p->xi*in_p.Cp/PI) + RG_MS[i]);
out_p->RZT_HC[i] = TZ[i];
out_p->GZ_HC[i] = GZ[i];
out_p->TT_HC[i] = OBS_TRANS[i];
*/


/*Commented By Ting for adjusting the output canopy reflectance, Intensity may be necessary for vertical crown shape.
out_p->ro[i]= out_p->PG*(in_p.RG[i]*(1-out_p->xi*in_p.Cp/PI) + RG_MS[i]) + out_p->PT*(in_p.RT[i]*Intensity + RT_MS[i]) + out_p->ZG*GZ[i] + out_p->ZT*TZ[i]
+ TT[i]*(out_p->QQ2B*(1-out_p->Pti) + out_p->QQ1B*out_p->Pti);
         //End of the comment*/

out_p->ro[i]= out_p->PG*(in_p.RG[i]+ RG_MS[i]) + out_p->PT*(in_p.RT[i]+ RT_MS[i]) + out_p->ZG*GZ[i] + out_p->ZT*TZ[i]; //Added By Ting


//out_p->M_factor[i] = (out_p->ro[i]-(out_p->RG_HC[i]*out_p->PG))/(out_p->FOLIAGE_REF[i]*out_p->PT);

/*A new M_factor defined by Ting*/
//out_p->M_factor2[i] = (out_p->ro[i]-(out_p->RG_HC[i]*out_p->PG)-(out_p->FOLIAGE_REF[i]*out_p->PT))/(out_p->FOLIAGE_REF[i]*out_p->ZT);
/*out_p->M_factor[i] = (out_p->ro[i]-(out_p->RZT_HC[i]*out_p->ZT))/(out_p->RT_HC[i]*out_p->PT);*/

}


} else /* no foliage (LAI=0) */
{


for (i=wave_index_min;i<wave_index_max;i++)
{
out_p->ro[i] = (1+fd[i])*in_p.RG[i];
}

}
}
