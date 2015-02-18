/*******************************************************************/
/*  pq_sub.c ->void PQ_SUB(in,out,shape)                           */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine that allows some height dependence for the sunlit   */
/*  crown proportion by separating the crown in cone and cylinder  */ 
/*  Latest update                                February 4, 1997  */
/*******************************************************************/

 
# include <math.h>
# include <stdio.h> 
# include "data.h" 


void PQ_SUB(in_p,out_p,shape)

struct PARAMETER in_p; 
struct RESULT *out_p ;
int shape;

{
 
	double F,AA,AAC,QV,QS,QVC,QSC,X;

	X = atan(2*in_p.R/out_p->E_r);
	if (X < 0) X = -X;
 
	F = 1-1.*out_p->phi/X ; 

	QV =1-out_p->Pv;
	QS =1-out_p->Ps;   
    
	QVC =1-out_p->Pvc; 
	QSC =1-out_p->Psc; 

	if (out_p->Pvc == 1) QVC = 0.00000000000001;

 /**************** this is on the principal plane only **************/
	if(QS<=QV)   AA =QS;
	if(QSC<=QVC)  AAC = QSC;
	if (QS>QV) AA =QV;
	if (QSC>QVC) AAC =QVC;

 /********************other planes *********************************/ 

	if (F <0) F=0;
	if (F >1) F=1;

           
	if (out_p->phi > X ) F = 0;

	switch(shape)
	{
		case 1:  out_p->Pti = ((QVC*QSC+(AAC-QVC*QSC)*F)*out_p->tic+
                   ((QV)*(QS)+(AA-(QV)*(QS))*F)*out_p->tib )
                        /(out_p->tac*QVC +QV*out_p->tab);
		break;
 
		case 2:  out_p->Psc =0;
				 out_p->Pvc =0;
				 out_p->Pti = (out_p->tib+out_p->tic)/(out_p->tab+out_p->tac);
		break;
 
	}



}
