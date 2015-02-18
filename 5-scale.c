/********************************************************************/
/*  5-Scale 1.5                                                     */
/*  4-Scale Theory :                                                */
/*  Jing M. Chen,                        chenj@geog.utoronto.ca     */
/*  Sylvain G. Leblanc          Sylvain.Leblanc@ccrs.nrcan.gc.ca    */
/*  4-Scale Code:													*/
/*  Sylvain G. Leblanc                                              */
/*  (C) CCRS 2004                                                   */
/********************************************************************/
/*  New loop mode									   October 2004 */
/********************************************************************/

 

/***********************  C librairies *******************************/
# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>

/*********************** 5-Scale Librairy ****************************/


# include "data.h"          /* contains all the relevant structures  */ 

/*****************************main program ***************************/ 
main (argc,argv)
int argc;
char *argv[];                                /* parameter input file */ 
{
typedef  double tableau[NN]; 
struct PARAMETER in_p; 
struct RESULT  out_p; 

/*********************** declarations of subroutine ********************/
 

double  exp();                   
double  pow();                         
double  log();
 
double  overlap();                     /* overlap function subroutine */
void Initialise();
void LIBERTY(int Default,struct PARAMETER in_p, struct RESULT *out_p); 
void FOUR_SCALE(struct PARAMETER in_p, struct RESULT *out_p);
void OPTICAL(struct PARAMETER in_p, struct RESULT *out_p);
void MULTIPLE_SCATTERING(struct PARAMETER in_p, struct RESULT *out_p);
void GetParameters(struct PARAMETER *in_p); 

/*******************  declaration of in/out pointers *******************/

FILE    *fp1,*fp2;                   /* for saving simulated data in a file */
char *com_file="longer_default_name";       /* input parameters file   */
int first_time=1;                  /* used to print name of parameters */


/************************** variables declaration *************************/
/**************************      integer          *************************/

int     i,k,kk,kkk,back,ref,trans;         /*  index of arrays  */
int		i_d,i_l,i_Omega,i_Gamma,i_m2,i_hb,i_ha,i_r,i_L,i_B,i_DOMAIN,i_SHAPE,i_Q; /* index of loop array */
int		file_index;
int     ii=0;  
int		NumberOfAngleRead=0;  
double  MU=0,V=0;
double OMEGA_T=0;
int i_wave=0;
int aaa=0;
int Dummy;

/**************************************************************************/
/******************* begining of the main program *************************/
/**************************************************************************/

printf("\n **************************************************** ");
printf("\n * 5-SCALE 1.5   (C) CCRS 2004 by  Chen and Leblanc *");
printf("\n *            Natural Resources Canada              *");
printf("\n * Problem or info: S. Leblanc 450-926-4646         *");
printf("\n **************************************************** \n");

/* first read some default values, some of which will be overwritten in loops */

Initialise(&in_p);

printf("\n Last modified October 8, 2004");


in_p.SPECTRAL=2; /* Multi spectral, from an old version, I am not sure it is still needed */


out_p.SZA_TMP=-1;  
out_p.DIST=1;                     /* Put this in the loop if you 
                                    need to change the tree distribution */ 


/* assign reflectivity to reflectivity array, only need to do this once */

i=0;
aaa=0;
printf("\n");



/* read parameter in input file in_p->COM_FILE, see Initialise.c, default:  "Input5Scale.txt"); */

GetParameters(&in_p);



printf("\n ****************************************************\n");

/* assigning optical inputs, this is an important part to change to run the model in hyperspectral mode */

for(i_wave =0; i_wave <2; i_wave++)
{
	for(back=0;back<in_p.NN_OPTIC;back++)	
	{
		for(ref=0;ref<in_p.NN_OPTIC;ref++)
		{
			for(trans=0;trans<in_p.NN_OPTIC;trans++)  
			{
				if(i_wave==0) /* Red */
				{
					out_p.BACKGROUND_REF[i]=in_p.OPTIC_REDG[back];
					out_p.FOLIAGE_REF[i]=in_p.OPTIC_REDT[ref];
					out_p.FOLIAGE_TRANS[i]=in_p.OPTIC_REDTT[trans];
					out_p.Wave[i] = 670;  /* used for Rayleight scattering in multiple scattering scheme */
					
					if((out_p.FOLIAGE_REF[i] + out_p.FOLIAGE_TRANS[i])<1) 
					{
						
						printf("Red: %f %f %f\t",out_p.BACKGROUND_REF[i], out_p.FOLIAGE_REF[i],out_p.FOLIAGE_TRANS[i]);
						i++;
						aaa++;
						if(aaa>1)
						{
							printf("\n");
							aaa=0;
						}
						
					}
				}
				else  /* NIR */
				{
					out_p.BACKGROUND_REF[i]=in_p.OPTIC_NIRG[back];
					out_p.FOLIAGE_REF[i]=in_p.OPTIC_NIRT[ref];
					out_p.FOLIAGE_TRANS[i]=in_p.OPTIC_NIRTT[trans];
					out_p.Wave[i] = 800; /* used for Rayleight scattering in multiple scattering scheme */
					
					if((out_p.FOLIAGE_REF[i] + out_p.FOLIAGE_TRANS[i])<1) 
					{
						printf("NIR: %f %f %f\t",out_p.BACKGROUND_REF[i], out_p.FOLIAGE_REF[i],out_p.FOLIAGE_TRANS[i]);
						i++;
						aaa++;
						if(aaa>1)
						{
							printf("\n");
							aaa=0;
						}					
					}
				}
			}
		}
	}
}
out_p.Optic_index = i;
printf("\n%d optical cases\n",out_p.Optic_index);


in_p.ALPHA_B=in_p.ALPHA_B*PI/180.;
in_p.ALPHA_L=in_p.ALPHA_L*PI/180.;
in_p.ALPHA = in_p.ALPHA*PI/180.;  /* January 12, 2004 */

file_index =0; /* this value can be changed to start output file numbering at otehr number than 0 */

Dummy = 0;

for(i_d=0;i_d<in_p.NN_D;i_d++)  
{
	in_p.D = in_p.D_RANGE[i_d];

	for(i_l=0;i_l<in_p.NN_LAI;i_l++)   
	{

		in_p.LAI = in_p.LAI_RANGE[i_l];


		for(i_Omega=0;i_Omega<in_p.NN_OMEGA;i_Omega++)
		{
			in_p.OMEGA_E = in_p.OMEGA_RANGE[i_Omega];


			for(i_Gamma=0;i_Gamma<in_p.NN_GAMMA;i_Gamma++)
			{

				in_p.GAMMA_E = in_p.GAMMA_RANGE[i_Gamma];


				for(i_m2=0;i_m2<in_p.NN_M2;i_m2++)
				{
					in_p.m2 = in_p.M2_RANGE[i_m2];

					for(i_hb=0;i_hb<in_p.NN_HB;i_hb++)
					{
						in_p.Hb = in_p.HB_RANGE[i_hb];


						for(i_ha=0;i_ha<in_p.NN_HA;i_ha++)
						{
							in_p.Ha = in_p.HA_RANGE[i_ha];


							for(i_r=0;i_r<in_p.NN_R;i_r++)
							{
								in_p.R = in_p.R_RANGE[i_r];
								out_p.DIST=1;

								for(i_L=0;i_L<in_p.NN_ANGLE;i_L++)
								{
									in_p.ALPHA_L = in_p.ALPHA_G[i_L];

									for(i_B=0;i_B<in_p.NN_ANGLE;i_B++)
									{
										in_p.ALPHA_B = in_p.ALPHA_G[i_B];

										if(in_p.ALPHA_B<0 && in_p.ALPHA_L<0) 
										{
												(in_p.GE_CHOICE,"NO_BRANCH"); 
												
										}
										else strcpy(in_p.GE_CHOICE,"BRANCH");

										for(i_DOMAIN=0;i_DOMAIN<in_p.NN_B;i_DOMAIN++)
										{
											in_p.B = in_p.B_RANGE[i_DOMAIN];

											for(i_SHAPE=0;i_SHAPE<in_p.NN_SHAPE;i_SHAPE++)
											{

												if(in_p.SHAPE_RANGE[i_SHAPE] == 1.)
												{
													strcpy(in_p.SHAPE,"CONE_CYLINDER"); 
												} 
												else if(in_p.SHAPE_RANGE[i_SHAPE] == 2.)
												{
													strcpy(in_p.SHAPE,"SPHEROID");
												} 
												else
												{
													printf("\n Problem with shape %f,(%d)", in_p.SHAPE_RANGE[i_SHAPE],i_SHAPE);
													exit(0);
												}
	
				
												for(i_Q=0;i_Q<in_p.NN_QUADRAT;i_Q++)
												{

													in_p.n = in_p.QUADRAT_RANGE[i_Q];

													

													if	(!strcmp(in_p.SHAPE,"CONE_CYLINDER")) V = PI*in_p.R*in_p.R*in_p.Hb;  /* cylinder approximation */
													else  V = 2/3.*PI*in_p.R*in_p.R*in_p.Hb;  

													MU =  in_p.LAI*in_p.B/(in_p.D*V);
													in_p.ALPHA_B=in_p.ALPHA_B*PI/180.;
													in_p.ALPHA_L=in_p.ALPHA_L*PI/180.;

													/* the next line creates the filename */
													sprintf(in_p.OUTPUT_FILE,"5SCALE_OUTPUT_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d.txt",
															file_index,i_d,i_l,i_Omega,i_Gamma,i_m2,i_hb,i_ha,i_r,i_L,i_B,i_DOMAIN,i_SHAPE,i_Q);

													printf("\n%s",in_p.OUTPUT_FILE);

													//  1) file_index,
													//  2) i_d,		TREE DENSITY
													//  3) i_l,		LAI
													//  4) i_Omega, Within crown Clumping
													//  5) i_Gamma, Needle-to-shoot 
													//  6) i_m2,	Neyman
													//  7) i_hb,	Height of crown
													//  8) i_ha,	"stick" height
													//  9) i_r,		Crown radius
													// 10) i_L,		Foliage orientation 
													// 11) i_B, 	Branch Orientation
													// 12) i_DOMAIN Domain size 
													// 13) i_SHAPE  cone+cylinder or spheroid
													// 14) i_Q		Nuber of quadrat (default 40)


												
													//if(MU > 0.001 && MU < 20 && (PI*in_p.R*in_p.R*in_p.D < 5.0*in_p.B)) // only simulate when crown are with mu = 0.001 and 20 and crown do not overlap more than 5 times
												
													if((PI*in_p.R*in_p.R*in_p.D < 5.0*in_p.B))
													{

														if(!(fp2=fopen(in_p.OUTPUT_FILE,"r")))
														{
														

															fp1=fopen(in_p.OUTPUT_FILE,"w");

															file_index++;

												
															for(kk=0;kk<in_p.NN_SZA;kk++)for(kkk=0;kkk<in_p.NN_PHI;kkk++)for (k=0;k<in_p.NN_VZA;k++) 
															{ 
																
  																out_p.vza = in_p.VZA[k]*PI/180.;

																out_p.phi = in_p.PHI[kkk]*PI/180.;
																in_p.SZA = in_p.SZA2[kk]*PI/180.;

														
																if(in_p.LAI>0)
																{

																	/***************************************************************************/
   																	FOUR_SCALE(in_p,&out_p);          /* call to main subroutine in 4-scale3.c */
																	/***************************************************************************/
															
																	OMEGA_T = log(out_p.Pvg)/(log(exp(-out_p.GFoliage*in_p.LAI/cos(out_p.vza))));

																	//printf("%6.3f",OMEGA_T);

																	fprintf(fp1,"%4.1f \t%4.1f \t%4.1f",out_p.vza*180./PI,out_p.phi*180/PI,in_p.SZA*180/PI);

																/* the next line puts the sunlit and shaded proportions in the file */
																	fprintf(fp1,"\t%6.4f \t%6.4f \t%6.4f \t%6.4f \t%6.4f \t%6.4f\t", 
																		out_p.PT,out_p.PG,out_p.ZT,out_p.ZG,out_p.Pvg, OMEGA_T); 

																	for(i=0;i<out_p.Optic_index;i++) fprintf(fp1,"%7.5f\t",out_p.ro[i]);
																	fprintf(fp1,"\n");
																	
																}
															}
															fclose(fp1); 
															
														} else
														{
															file_index++;
															fclose(fp2);
															printf("\n%s already exist",in_p.OUTPUT_FILE);
														}
														
													

													} else printf("\n*xxxxxxxxxx**%6.3f %s  %6.3f",MU,in_p.OUTPUT_FILE,PI*in_p.R*in_p.R*in_p.D/in_p.B);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}
          	
			
 

}  /* end of all */
 
 
