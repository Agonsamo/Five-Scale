/*******************************************************************/
/*  parameters.c ->void GetParameters(in)                          */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that reads the input parameter file                 */
/*  Latest update                                October 7, 2004   */
/*******************************************************************/



# include <stdio.h> 
# include <stdlib.h>
# include <string.h>
# include "data.h"

void GetParameters(struct PARAMETER *local_in_p)

{

FILE *fp; 
char in_str[20];
double para_str;
int i;

printf("\n");

if (!(fp=fopen(local_in_p->COM_FILE,"r")))
{ 
	fclose(fp);	

} else while(!feof(fp))

{
	
	
		//old input file

		fscanf(fp,"%s ",in_str);
		fscanf(fp,"%lf\n",&para_str);

		printf("%s   \t%5.0f\n", in_str, para_str);

		if(!strcmp(in_str,"N_VZA"))
		{
			local_in_p->NN_VZA = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->VZA[i]);
		}
		else if(!strcmp(in_str,"N_SZA"))
		{
			local_in_p->NN_SZA = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->SZA2[i]);
		}
		else if(!strcmp(in_str,"N_PHI"))
		{
			local_in_p->NN_PHI = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->PHI[i]);

		}
		else if(!strcmp(in_str,"N_LAI"))
		{
			local_in_p->NN_LAI=para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->LAI_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_OMEGA"))  
		{
			local_in_p->NN_OMEGA = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OMEGA_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_D")) 
		{
			local_in_p->NN_D = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->D_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_B")) 
		{
			local_in_p->NN_B = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->B_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_QUADRAT")) 
		{
			local_in_p->NN_QUADRAT = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->QUADRAT_RANGE[i]);

		}

		else if(!strcmp(in_str,"N_GAMMA"))
		{
			local_in_p->NN_GAMMA = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->GAMMA_RANGE[i]);

		}

		else if(!strcmp(in_str,"N_M2"))
		{
			local_in_p->NN_M2 = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->M2_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_HB"))	
		{
			local_in_p->NN_HB = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->HB_RANGE[i]);

		}

		else if(!strcmp(in_str,"N_HA"))
		{
			local_in_p->NN_HA = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->HA_RANGE[i]);

		}
		else if(!strcmp(in_str,"N_R"))
		{
			local_in_p->NN_R = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->R_RANGE[i]);
				
		}
		else if(!strcmp(in_str,"N_SHAPE"))
		{
			local_in_p->NN_SHAPE = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->SHAPE_RANGE[i]);
		

		}
		else if(!strcmp(in_str,"N_ANGLE"))	
		{
			local_in_p->NN_ANGLE = para_str;
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->ALPHA_G[i]);

		}
	
		else if(!strcmp(in_str,"N_OPTIC"))
		{
			local_in_p->NN_OPTIC = para_str;
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_NIRG[i]);
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_NIRT[i]);
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_NIRTT[i]);
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_REDG[i]);
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_REDT[i]);
			fscanf(fp,"%*s");
			for(i=0;i<para_str;i++)fscanf(fp,"%lf",&local_in_p->OPTIC_REDTT[i]);


		}

	}


   fclose(fp);
}


