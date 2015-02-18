/*******************************************************************/
/*  tree_size.c ->void TREE_SIZE(in,out,choice)                    */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that reduces the tree size when the quadrat density */
/*  is larger than the mean density                                */
/*  Latest update                                February 7, 1997  */
/*******************************************************************/

# include <math.h>
# include <string.h> 
# include <stdio.h>
# include "data.h"

void TREE_SIZE(in_p,out_p,CHOICE) 


struct PARAMETER in_p;
struct RESULT *out_p;
char *CHOICE;

{

	int i=0;
	double D=0,Vg=0,Vg_0=0;

	D= in_p.D/in_p.n; 
	if(!strcmp(CHOICE,"VZA"))
	{
		Vg_0=out_p->Vg_0; 
	}
	if(!strcmp(CHOICE,"SZA"))
	{
		Vg_0=out_p->Sg_0; 
	}

	for (i=0;i<NN;i++) 
	{
		if (i > D)
		{
			Vg = Vg + out_p->Px[i]*Vg_0*D/i;
		} else
		{
			Vg = Vg + out_p->Px[i]*Vg_0;
		}   
	} 
	if(!strcmp(CHOICE,"VZA"))
	{
		out_p->Vg=Vg;
	}
	if(!strcmp(CHOICE,"SZA"))
	{
		out_p->Sg=Vg;
	}
}  
