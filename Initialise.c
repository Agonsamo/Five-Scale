/*******************************************************************/
/*  Initialise.c ->void Initialise(in,out)                         */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  This subroutine  initialise the input parameter                */
/*  Latest update                                  March 6, 2000   */
/*******************************************************************/

# include <stdio.h> 
# include <stdlib.h>
# include <string.h>
# include "data.h"

void Initialise(struct PARAMETER *in_p)
{


	strcpy(in_p->ANGLE_FILE,"angle.txt");  /* old default, not used anymore */
	strcpy(in_p->OUTPUT_FILE,"output.txt");
	strcpy(in_p->COM_FILE,"Input5Scale.txt");
/*	strcpy(in_p->OPTICAL_FILE,"optical.txt"); */
	in_p->SPECTRAL   =	1;		
	strcpy(in_p->GE_CHOICE,"NO_BRANCH");	
/*	strcpy(in_p->SHAPE,"CONE_CYLINDER"); */
	strcpy(in_p->SHAPE,"SPHEROID");
	in_p->Ha   =  1.0;	
	in_p->Hb   =  5.0;	
	in_p->A    =  0.00;	
	in_p->C    =  0.50;	
	in_p->Cp   =  1.0;
	in_p->LAI  =  3.50;	
	in_p->B    =  10000.0;	
	in_p->D    =  1000;
	in_p->n    =  40;	
	in_p->R    =  1.0;	
	in_p->m2   = 2;	
	in_p->SZA  = 45.0; 	
	in_p->Band[1] = 670.0; /* POLDER red band */
	in_p->Band[2] = 865.0; /* POLDER NIR band */
	in_p->Band[3] = 1600.0; /* SWIR */
	in_p->Band[4] = 1600.0; /* SWIR */
	in_p->G1   =  0.050;
	in_p->GZ1  =  0.001;	
	in_p->G2   =  0.270;
	in_p->GZ2  =  0.010;
	in_p->G3   =  0.200;
	in_p->GZ3  =  0.005;
	in_p->G4   =  0.200;
	in_p->GZ4  =  0.005;
	in_p->T1   =  0.070;
	in_p->TZ1  =  0.001;
	in_p->T2   =  0.470;
	in_p->TZ2  =  0.010;
	in_p->T3   =  0.100;
	in_p->TZ3  =  0.005;
	in_p->T4   =  0.100;
	in_p->TZ4  =  0.005;
	in_p->TT1  =  0.020;
	in_p->TT2  =  0.300;
	in_p->TT3  =  0.150;
	in_p->TT4  =  0.150;
	in_p->Ws   =  0.40 ;			
	in_p->OMEGA_E  = 0.98;			
	in_p->GAMMA_E  = 1.41;
	in_p->ALPHA_B  = 10.0;			
	in_p->ALPHA_L  = 20.0;			
	in_p->Ll       = 0.80;				
	in_p->Fr       = 0.00;			
	in_p->ALPHA    = 13.0;			
	in_p->RATIO    = 0.20;			
	in_p->Rb       = 0.10;			
	in_p->DeltaLAI = 0.20;	

	in_p->LIBERTY_DEFAULT =	1;

	strcpy(in_p->PIGMENT_FILE,"pigment.txt");
	strcpy(in_p->WATER_FILE,"water.txt"); 
	strcpy(in_p->ALBINO_FILE,"albino.txt");
	strcpy(in_p->LIGCELL_FILE,"ligcell.txt");
	strcpy(in_p->PROTEIN_FILE,"protein.txt");
		
	in_p->m_D        = 40.0;			
	in_p->m_XU       = 0.045;			
	in_p->m_THICK    = 1.600;			
	in_p->m_BASELINE = 0.0005;			
	in_p->m_ELEMENT  = 2.000;			
	in_p->m_C_FACTOR = 200.0;			
	in_p->m_L_FACTOR = 40.00;			
	in_p->m_P_FACTOR = 1.000;			
	in_p->m_W_FACTOR = 100.0;		


}