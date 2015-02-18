/*******************************************************************/
/*  xi.c ->double XI(sza,vza,phi)                                  */
/*  Jing M. Chen,                     jing.chen@ccrs.nrcan.gc.ca   */
/*  Sylvain G. Leblanc          sylvain.leblanc@ccrs.nrcan.gc.ca   */
/*******************************************************************/
/*  Subroutine thatcalculte the scattering angle (phase angle)     */
/*  Latest update                                November 30,1995  */
/*******************************************************************/


# include <math.h>
# include <stdio.h>


double XI(sza,vza,phi)
double sza,vza,phi;

{

  double Uy,Uz;
  double Vy,Vz;
  double answer;
 

  Uy = sin(sza);
  Uz = cos(sza);


  Vy = sin(vza)*cos(phi);
  Vz = cos(vza);


  answer = (Uy*Vy+Uz*Vz);

  answer = acos(answer);
  return answer ;

}
