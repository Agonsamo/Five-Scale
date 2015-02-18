/*************************************************************************/
/*  4-scale  (version  3.3 for 5-Scale 5.2)      void FOUR_SCALE()       */
/*************************************************************************/
/*  Jing M. Chen, @CCRS                       jing.chen@ccrs.nrcan.gc.ca */
/*  Sylvain G. Leblanc @CCRS            sylvain.leblanc@ccrs.nrcan.gc.ca */
/* (C) 2000 Canada Centre for Remote Sensing	                         */
/*  588 Booth Street, 4th Floor, Ottawa, Ontario, Canada, K1A 0Y7        */
/*************************************************************************/

/*********************** usual C librairies ******************************/

# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include <string.h>

/*********************** 4-Scale header file *****************************/ 

# include "data.h" 


void FOUR_SCALE( in_p, out_p)

struct PARAMETER in_p;
struct RESULT *out_p;		
{
 
	/* Standard C function declaration   */
	double  exp();                     	
	double  pow();         
	double  log();
	double  cos();
	double  sqrt();
	double  asin();
 
	/* declaration of functions used by FOUR_SCALE() */
	double  XI();                           /* phase or scattering angle */
	void    DISTANCE();
	void    FO();
	void    LS();
	void    MULTIPLE_SCATTERING();
	void    OVERLAP();
	void    P_GAP_AX_B();
	void    P_GAP_BRANCH();
	void    PTG_SUB();
	void    PQ_SUB();
	void    PS();
	void	POISSON();	               /* random (Poisson) distribution  */
	void    Q();
	void    NADIR();
	void	NEYMAN();              /* clumped distribution using Neyman  */
	void	CONE_TA();
	void    CONE_TI();
	void	CONE_GS();
	void	SPHEROID_TA();
	void    SPHEROID_TI();
	void    SPHEROID_SVG();
	void    TREE_SIZE();

 
	/********************* variables declaration *************************/
	/*********************      integer          *************************/
	int     k=0;                                    /*  index of arrays  */
	int     distribution =2;               /* 1 = Poisson ,  2 = neyman  */
	int     Gchoice=1;                   /* 1--> G()=a+c(), 2-> branches */
	int     shape=1;                /* 1--> cone+Cylinder 2--> spheroid  */


	/*********************      double           ************************/


	double	temp_Pig=0;
	double  temp_Px[NN];
	double	num=0,dem=0;
	double  Intensity=0;
	double  TT1=0,TT2=0,TT3=0,TT4=0;





	out_p->OmegaE_stand = in_p.OMEGA_E;
 
 

	/********* calculating the height of the cone part for conifers ******/
 
	if(tan(in_p.ALPHA)!=0)out_p->Hc=in_p.R/tan(in_p.ALPHA);
	else out_p->Hc =0;

	
	if (!strcmp(in_p.SHAPE,"SPHEROID"))
	{
		out_p->Hc =0;      /* no cone for the spheroid */
		shape=2;
	}
	else shape=1;



	if (!strcmp(in_p.GE_CHOICE,"BRANCH")) Gchoice=2; 
	else Gchoice=1;



	/****  Neyman grouping m2=0 forces a Poisson distribution    ****/ 

	if(in_p.m2==0) distribution =1;   /* Poisson tree distribution  */
	if(in_p.m2 >0) distribution =2;   /* Neyman tree distribution   */
 

	/********* call to subroutine of trees distribution *************/
  
     
	out_p->Error[4]=0; /* error initialisation, NT/95 only  */
	out_p->Error[1]=0; /* error initialisation, NT/95 only  */

	/* should change the next lines so that out_p.Error[] can be use with 
	Command line version (UNIX/DOS) */

	if(PI*in_p.R*in_p.R*in_p.D > in_p.B)
	{
		if(in_p.Fr >0.5 && out_p->DIST==1)
		{
			printf("\n Maybe you should reduce the repulsion factor \n");
			printf("\n because the total surface of the crowns is larger than the domain surface");
			out_p->Error[4]=1;
		}
 		
		/* Need to introduce something 
		   to take care automatically of 
		   this repulsion problem 
		   (done in Windows 95/NT version) */
	}


	/* calculation of scattering angle */

	out_p->xi = XI(in_p.SZA,out_p->vza,out_p->phi);
 
 

  /********************************************************************/
  /********* call to subroutines for shaded and sunlit areas **********/
  /********************************************************************/
 
 
	switch (shape)
	{
 
		case 1: CONE_TA(in_p,out_p);
			    CONE_GS(in_p,out_p,"VZA");
                if(out_p->SZA_TMP!=in_p.SZA || (out_p->DIST==1))
				{		  
					CONE_GS(in_p,out_p,"LAI");
					CONE_GS(in_p,out_p,"SZA");
				}
		break;
 
		case 2: SPHEROID_TA(in_p,out_p);
				SPHEROID_SVG(in_p,out_p,"VZA");
				out_p->Vgc=0;
				if(out_p->SZA_TMP!=in_p.SZA || (out_p->DIST==1))
				{
			
					SPHEROID_SVG(in_p,out_p,"LAI");
					SPHEROID_SVG(in_p,out_p,"SZA");
					out_p->tic=0;
					out_p->Sgc=0;
				}
				out_p->tac=0.;
		break;
 
		default: printf("\n error, bad tree shape requested");
        exit(1);
 
	}
 
 
	if (out_p->xi <= 0.000001 && out_p->xi >=0.0)
	{
		out_p->tic = out_p->tac;
		out_p->tib = out_p->tab;
	}
	else
	{
		switch (shape)
		{
			 case 1: CONE_TI(in_p,out_p);
 
			 break;
 
			 case 2: SPHEROID_TI(in_p,out_p);
 
			 break;
 
			 default: printf("\n error, bad tree shape requested");
			 exit(1);
 
		 }
	}
 

	out_p->Sv = out_p->V/(out_p->Vg_0*cos(out_p->vza));
	out_p->Ss = out_p->V/(out_p->Sg_0*cos(in_p.SZA));


	/****************************************************************/
    /*  out_p->DIST=1; <--USE ONLY INSIDE A LOOP THAT CHANGES       */
    /*  EITHER THE TREE DENSITY,QUADRAT SIZE, OR NEYMAN GROUPING    */
    /*  For Plane plots, the program only compute Px once           */
    /****************************************************************/

	if(out_p->DIST==1)
	{
		for (k=0;k<NN;k++) out_p->Px[k] =0.;
		switch(distribution)
		{
			case 1: POISSON(in_p,out_p);
			break;
 
			case 2: NEYMAN(in_p,out_p);
			break;
 
			default: printf("\n error, bad tree distribution requested");
			exit(1);
		 }
    
	} 	  

	TREE_SIZE(in_p,out_p,"VZA") ;

	if(out_p->Pvg< 0) out_p->Pvg =0;

	switch (Gchoice)
	{
		case 1:  P_GAP_AX_B(in_p,out_p,"VZA");
            
				 out_p->PSG0_VIEW = 1- out_p->Vg*in_p.D/in_p.B*(1-out_p->PgapV);
 
				 if(out_p->SZA_TMP!=in_p.SZA)
				 {	  
					 P_GAP_AX_B(in_p,out_p,"LAI");
					 P_GAP_AX_B(in_p,out_p,"SZA");
					 TREE_SIZE(in_p,out_p,"SZA") ;
					 out_p->PSG0_SUN = 1- out_p->Sg*in_p.D/in_p.B*(1-out_p->PgapS);
					 P_GAP_AX_B(in_p,out_p,"0");
				 }
		break;
 
		case 2:
				 
				 P_GAP_BRANCH(in_p,out_p,"VZA");
				 out_p->PSG0_VIEW = 1- out_p->Vg*in_p.D/in_p.B*(1-out_p->PgapV);
 
				 if(out_p->SZA_TMP!=in_p.SZA)
				 {
					 P_GAP_BRANCH(in_p,out_p,"LAI");
					 P_GAP_BRANCH(in_p,out_p,"SZA");
					 TREE_SIZE(in_p,out_p,"SZA") ;
					 out_p->PSG0_SUN = 1- out_p->Sg*in_p.D/in_p.B*(1-out_p->PgapS);
					 P_GAP_AX_B(in_p,out_p,"0");
				 }
		break;
	}
 
	out_p->PSG_HOT0= 1- PI*in_p.R*in_p.R*in_p.D/in_p.B*(1-out_p->Pgap0);
 


	if (out_p->PSG0_SUN  <0) 
	{
		/* when this condition is met, it means that the projection*number of trees is larger than the domain */ 
		out_p->PSG0_SUN =0;	
	}
	if (out_p->PSG0_VIEW <0) 
	{
		/* when this condition is met, it means that the projection*number of trees is larger than the domain */ 
		out_p->PSG0_VIEW =0;
	}


	 
	if (out_p->SZA_TMP != in_p.SZA || (out_p->DIST==1) )
	{
		
		for (k=0;k<NN;k++) 
		{
			temp_Px[k] = out_p->Px[k] ;
			out_p->Px[k] =0;
		}

		POISSON(in_p,out_p);       /* call POISSON for computation of Omega_T  ... */
		temp_Pig = out_p->Pig;
		out_p->Pig=0;   
		OVERLAP(in_p,out_p,"SZA"); /* Pig (Poisson) Eq. 43 [1] */

		out_p->Pig_poisson = out_p->Pig;
		if(out_p->Pig_poisson >1 || out_p->Pig_poisson<0) out_p->Pig_poisson=0;
		out_p->Pig = temp_Pig;
		for (k=0;k<NN;k++) 
		{
			out_p->Px[k] =temp_Px[k] ;
		}
		
	}


	/*****************************  Pig  ********************************/
	/********* probability of having an illuminated bacground ***********/
	/********************************************************************/
 

	if (out_p->SZA_TMP != in_p.SZA || (out_p->DIST==1))
	{
 

		OVERLAP(in_p,out_p,"LAI");

		out_p->Pig=0;
		OVERLAP(in_p,out_p,"SZA");
		OVERLAP(in_p,out_p,"NADIR");


		FO(in_p,out_p,"SZA");


	/*****************************************************************/
	/******************* Calculation of OmegaT   *********************/
	/*****************************************************************/

		out_p->OmegaT = log(out_p->Pig)/log(out_p->Pig_poisson);
	}

	/*******************************  Pvg  ***************************/
	/************** probability of viewing the background ************/
	/*****************************************************************/
 
 
	
 
	out_p->Pvg=0;
	OVERLAP(in_p,out_p,"VZA");
	out_p->Error[13] =0;
	FO(in_p,out_p,"VZA");
	if(out_p->Pvg < 0 )
	{
		out_p->Pvg =0;
		out_p->Error[13] =1;
	}
 		  


	if(out_p->DIST==1) 
	{
		DISTANCE(in_p,out_p);
 

	}
 
 
	PS(in_p,out_p);    /* computes PS and parameters for PG */


	/************************* overlapping subroutine *******************/
 
	PQ_SUB(in_p,out_p,shape);


	out_p->Error[17] =0;

	if(out_p->Pti>1)
	{
		out_p->Pti=1;
		out_p->Error[17] =1;
	}
	else if(out_p->Pti<0) 
	{
		out_p->Pti=0;
		out_p->Error[17] =1;
	}

	/**********************************************************************/
	/*       calculation of PG --- Prob seen sunlit background            */
	/**********************************************************************/

	if(out_p->DIST==1)NADIR(in_p,out_p);

	
	PTG_SUB(in_p,out_p,"NADIR",120000,0.10); 
 
	PTG_SUB(in_p,out_p,"GROUND",20000,0.01) ;
 
 
 
 
	/**************************************************************************/
	/*              calculation of PT --- Prob seen sunlit foliage            */
	/**************************************************************************/

	num = log(out_p->Pvg_mean);

	dem = log(exp(-0.5*in_p.LAI/(0.537+0.025*in_p.LAI))); /* based on Chen et al, 1999 */

	out_p->OmegaTotal= num/dem;
 

 
	LS(in_p,out_p,shape);



	Q(in_p,out_p,in_p.Cp,shape);  


	/******************************  PT   *********************************/
 
	out_p->PT_Cold =  out_p->QQ1*out_p->Pti +(1-out_p->Pti)*out_p->QQ2  ;
 
	PTG_SUB(in_p,out_p,"CANOPY",40000,0.001) ;

	/***************** reflectance computation ****************************/
 
	out_p->Error[2] =0;
	out_p->ZG = out_p->Pvg - out_p->PG;


	if(out_p->PT<0) 
	{
		out_p->PT =0;
	}
	if (out_p->ZG <0)
	{
		out_p->ZG =0;
		out_p->PG = out_p->Pvg;
		if (out_p->xi > 0.0001) /* it doesn't matter if this is not calculated properly at the hotspot */
		{
		/*	printf("\n BRDF may not be calculated correctly");
			printf("(PG at vza = %5.1f deg.)\n",out_p->vza*180./PI);*/
			out_p->Error[2] =1;

		}
	}

	out_p->ZT = (1-out_p->Pvg) - out_p->PT;

	out_p->Error[3] =0;

	if (out_p->ZT <0)
	{
		out_p->ZT =0;
		out_p->PT = 1-out_p->Pvg;
		if(out_p->xi >0.0001)
		{
			/*printf("\n  BRDF may not be calculated correctly");
			printf("(PT at vza  =%5.1f deg.)\n",out_p->vza*180./PI);*/
			out_p->Error[3] =1;
		}
	}
 


	/********************  MULTIPLE SCATTERING  ****************/


	MULTIPLE_SCATTERING(in_p,out_p);
	

	
	out_p->SZA_TMP = -10000;
 
	out_p->DIST=0; 
	/* out_p->DIST=0; means that if DIST is not changed elsewhere, the tree distribution */
	/*  and other parameters won't be re-calculated for the next data point              */
	
}
 

