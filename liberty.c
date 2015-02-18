/*******************************************************************/
/*  liberty.c ->void LIBERTY(Default,in,out)                       */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine LIBERTY by T. Dawson                                */
/*  Latest update                                 March, 2000      */
/*******************************************************************/


# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include "data.h"
# include "liberty.h"



void LIBERTY(Default , in_p, out_p)

struct PARAMETER in_p;
struct RESULT *out_p;
int Default;

{

	
	FILE *infile;
	int i;
	int MAX_was_420 = 500;
	para_r=0;
	vert_r=0;
	
		D = in_p.m_D;
		xu = in_p.m_XU;
		thick = in_p.m_THICK;
		baseline = in_p.m_BASELINE;
		element = in_p.m_ELEMENT;
		c_factor = in_p.m_C_FACTOR;
		w_factor = in_p.m_W_FACTOR;
		l_factor = in_p.m_L_FACTOR;
		p_factor = in_p.m_P_FACTOR;
	
	



	/* reading in pigment absorption coefficients... */

	if(Default) 
	{
		for(i=0;i<MAX_was_420;i++)
		{
			ke[i] = 0;
			k_chloro[i] = 0;
			k_water[i]= WATER_DEFAULT[i];
		 
			k_ligcell[i] = LIGCELL_DEFAULT[i];
			k_protein[i] = PROTEIN_DEFAULT[i];
		}
		for(i=0;i<120;i++)
		{
			k_chloro[i] = PIGMENT_DEFAULT[i];
			ke[i] = ALBINO_DEFAULT[i];
		}

		


	} else
	{
		if ((infile = fopen(in_p.PIGMENT_FILE,"rt")) == NULL) 
		{
			printf("\nCould not Open File %s",in_p.PIGMENT_FILE);

		} else
		{
			rewind(infile);
			printf("Reading pigment absorption coefficients...");
			i = 0;
			while (!feof(infile)) 
			{
				 fscanf(infile,"%lf",&k_chloro[i]);
				 i=i+1;
			}
			fclose(infile);


			/* Ok so far, Now to read in the absorption data... */

			/* Now to set chloro absorption to zero for higher wavelengths */
			for (i = 120; i <= MAX_was_420; i++) 
			{
				k_chloro[i] = 0;
			}

			/* similarly for the water absorption coefficients... */

			if ((infile = fopen(in_p.WATER_FILE,"rt")) == NULL) 
			{
				printf("\nCould not Open File %s",in_p.WATER_FILE);
				exit(0);
			} else
			{
				/* Ok so far, Now to read in the absorption data... */
				rewind(infile);
				printf("Reading water absorption coefficients...");
				i = 0;
				while (!feof(infile)) 
				{
					 fscanf(infile,"%lf",&k_water[i]);
					 i=i+1;
				}
				fclose(infile);



				/* now to read in trace element absorption coefficients... */

				if ((infile = fopen(in_p.ALBINO_FILE,"rt")) == NULL) 
				{
					printf("\nCould not Open File water.dat");
					exit(0);
		
				} else
				{
					/* reset ke[i] array... */
					for (i = 0; i <= MAX_was_420; i++) 
					{
						ke[i] = 0;
					}
					/* Ok so far, Now to read in the absorption data... */
					rewind(infile);
					printf("Reading trace element absorption coefficients...");
					i = 0;
					while (!feof(infile)) 
					{
						 fscanf(infile,"%lf",&ke[i]);
						 /* printf("albino: %d %f\n",i,ke[i]); */
						 i=i+1;
					}
					fclose(infile);

					/* now to read in cellulose and lignin absorption coefficients... */

					if ((infile = fopen(in_p.LIGCELL_FILE,"rt")) == NULL) 
					{
						printf("\nCould not Open File water.dat");
						exit(0);
					}
					else
					{
						/* Ok so far, Now to read in the absorption data... */
						rewind(infile);
						/* first reset k_ligcell[i] array... */
						for (i = 0; i <= MAX_was_420; i++) 
						{ 
							k_ligcell[i] = 0; 
						}
						printf("Reading cellulose and lignin absorption coefficients...");
						i = 0; /* Note: this is because the data starts from 400 nm */
						while (!feof(infile)) 
						{
							 fscanf(infile,"%lf",&k_ligcell[i]);
							 /* printf("%d %f\n",(400+(i*5)),k_ligcell[i]); */
							 i=i+1;
						}
						fclose(infile);
					


						/* now to read in protein absorption coefficients... */
						if ((infile = fopen(in_p.PROTEIN_FILE,"rt")) == NULL) 
						{
							printf("\nCould not Open File water.dat");
							exit(0);
					
						} else
						{
							/* Ok so far, Now to read in the absorption data... */
							rewind(infile);
							printf("Reading protein absorption coefficients...");
							i = 0; /* Note: this is because the data starts from 400 nm */
							while (!feof(infile)) 
							{
								 fscanf(infile,"%lf",&k_protein[i]);
								 i=i+1;
							}
							fclose(infile);
						}
					}
				}
			}
		}
	
	}


	

	/* set the following loop for MAX_was_420 (to 2500nm) or 120 (1000nm) */
	for (i = 0; i <= MAX_was_420; i++) 
	{
		/* determining the total absorption coefficient... */
		coeff = (D * (baseline+(k_chloro[i]*c_factor)+(k_water[i]*w_factor)+(ke[i]*element)+(k_ligcell[i]*l_factor)+(k_protein[i]*p_factor)));
		/* change of refractive index over wavelength... */

		N1=1.4891-(0.0005*i);
		refrac();
		para_rad();
		vert_rad();
		eval_me();
		eval_mi();

		calc_M();
		calc_T();
		eval_x();
		calc_R();

		/* The next bit works out transmittance based upon Benford...*/
		/* setting up unchanging parameters... */
		a=(2*x*me)+(x*T)-(x*T*2*x*me);
		rb=a;
		tb=sqrt(((R-rb)*(1-(R*rb)))/R);
		contfunc();
	
	

		wavelength[i] = 400+(i*5);
		Refl_max[i]  = R;
		out_p->FOLIAGE_REF[i] = refl;
		out_p->FOLIAGE_TRANS[i]= trans;

	
	
	
			
	}
	


}

void refrac()
{
	/* average angle of incident light */
	in_angle = 59;
	/* Index of refraction */
	N0 = 1.0;
	alpha = in_angle * PI/180;
	beta = asin((N0/N1) * sin(alpha));
	/* printf("Alpha: %f Beta: %f\n",alpha,beta); */
}



// working out the horizontal (parallel to plane) component of reflected 
// radiation (Snells Law of Refraction) 

void para_rad()
{
	/* printf("Testing coeficients of reflection...\n"); */
	para_r = (tan(alpha - beta))/(tan(alpha + beta));
	/* printf("Parallel: %f",para_r); */
}


void vert_rad()
{
	vert_r = -(sin(alpha-beta))/(sin(alpha+beta));
	/* printf(" Vertical: %f\n",vert_r); */
}

/*--------------------------------------------------------------------*/

/* working out the total reflected radiation...*/
void tot_ref()
{
	double plus, dif;

	plus = alpha + beta;
	dif = alpha - beta;

	refl =  0.5 * ( ((sin(dif)*sin(dif))/(sin(plus)*sin(plus))) +
		((tan(dif)*tan(dif))/(tan(plus)*tan(plus))) );
}


/*-------------------------------------------------------------------------*/

/* The evaluation of the regular reflectance for diffuse incident */
/* radiation me for angles of alpha between 0 and PI/2  */
void eval_me()
{
	int a;
	double width;

	me = 0;
	width = PI/180;

	for (a = 1; a <= 90; a++) 
	{
		alpha = a * PI/180;
		beta = asin(N0/N1 * sin(alpha));
		tot_ref();
		me = me + (refl * sin(alpha) * cos (alpha) * width);
	}
	me=me*2;
   /* printf("me : %f\n",me); */
}

/*-------------------------------------------------------------------------*/

void eval_mi()
{
	int a;
	double mint,width;

	mi = 0;
	mint = 0;
	width = PI/180;
	critical = asin(N0/N1)*180/PI;
	/* printf("Critical angle: %f\n",critical); */
	for (a = 1; a <= (critical); a++) 
	{
		alpha = a * PI/180;
		beta = asin((N0/N1) * sin(alpha));
		tot_ref();
		mint = mint + (refl * sin(alpha) * cos (alpha) * width);
	}
	mi = (1 - (sin(critical*PI/180))*(sin(critical*PI/180))) + (2*mint);
   /* printf("mi : %f\n",mi); */
}


/*-----------------------------------------------------------------------*/



/*-----------------------------------------------------------------------*/

	/* Evaluation of M, the total radiation reaching
	the surface after one pass through the sphere */

void calc_M()
{
	  M = (2/(coeff * coeff))*(1-(coeff+1)*exp(-coeff));
}

/*-----------------------------------------------------------------------*/

void calc_T()
{
	  T = ((1-mi) * M)/(1-(mi * M));
}

/*-----------------------------------------------------------------------*/

void calc_R()
{
	double a,b,c,next_R;
	int iterations;

	a = (me * T) + (x * T) - me - T - (x * me * T);
	b = 1 + (x * me * T) - (2 * x * x * me * me * T);
	c = (2 * me * x * x * T) - (x * T) - (2 * x * me);

	R = 0.5;   /* initial guess */
	for (iterations = 1; iterations < 50; iterations++) 
	{
		/* simple iterative method... */
		next_R = -(a*(R*R)+c)/b;
		/* printf("alt_r: %f\n", next_R); */
		R = next_R;
	}


}

/*-----------------------------------------------------------------------*/

void eval_x()
{
	x = xu / (1 - (1 - (2*xu)) * T);
}

/*-----------------------------------------------------------------------*/

/* subroutine to determine the value of R and T as a continuous function
   of thickness - Benford, F., 1946, 'Radiation in a diffusing medium',
   Journal of the Optical Society of America, 36, 524-537. */
void contfunc()
{
	double fraction,top,bot1,bot2,cur_t,cur_r,prev_t,prev_r;
	int step, whole;

	/* little trick to seperate the fractional part from the real number... */
	whole=thick;
	fraction=thick-whole;

	/* The next bit works out the fractional value */
	/* for the interval between 1 and 2... */
	top=pow(tb,1+fraction)*(pow((pow((1+tb),2)-pow(rb,2)),(1-fraction)));
	bot1=(pow((1+tb),(2*(1-fraction)))-pow(rb,2));
	bot2=(1+((64/3)*fraction)*(fraction-0.5)*(fraction-1)*0.001);
	tif=top/(bot1*bot2);
	rif=(1+pow(rb,2)-pow(tb,2)-sqrt(pow((1+pow(rb,2)-pow(tb,2)),2)-(4*pow(rb,2)*(1-pow(tif,2)))))/(2*rb);

	/* Now to work out for integral thickness greater than 2 ... */
	if (whole >= 2) 
	{
		prev_t=1;
		prev_r=0;
		for (step = 1; step<= (whole-1); step++) 
		{
		   cur_t=(prev_t * tb)/(1-(prev_r*rb));
		   cur_r=prev_r + (((prev_t*prev_t)*rb)/(1-(prev_r*rb)));
		   prev_t=cur_t;
		   prev_r=cur_r;
		}
	}
	else 
	{
		cur_t=1;
		cur_r=0;
	}
	/* Combine the two results for thickness from 1 to infinity... */
	trans=(cur_t*tif)/(1-(rif*cur_r));
	refl=cur_r+((cur_t*cur_t*rif)/(1-(rif*cur_r)));
}


