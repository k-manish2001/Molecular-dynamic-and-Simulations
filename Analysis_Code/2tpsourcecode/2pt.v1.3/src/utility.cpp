/*********************************************************************************
*                     Two-Phase Thermodynamics (2PT) Program                     *
*                          Shiang-Tai Lin (stlin@ntw.edu.tw)                     *
* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *
*                     Prabal K. Maiti (maiti@physics.iisc.ernet.in)              *
*  Department of Physics, Indian Institute of Science, Bangalore, India, 560012  *
*                          Tod A Pascal (tpascal@wag.caltech.edu)                *
*                                     and                                        *
*                        William A Goddard III (wag@wag.caltech.edu)             *
*         Materials and Process Simulation Center, Caltech, Pasadena, CA USA     *
*                                Copyright (c) 2010                              *
*                                All rights Reserved                             *
*                             *************************                          *
*                                      Cite:                                     *
* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *
* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *
* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13(1), 169-181             *
***********************************************************************************/

/* General routines, some from Numerical Recipes (slightly modified). */
#include "utility.h"
#include <math.h>

void ROTATE(double **a,int i,int j,int k,int l,double *tau,double *s)
{
  double g,h;
  g=a[i][j];
  h=a[k][l];
  a[i][j]=g-(*s)*(h+g*(*tau));
  a[k][l]=h+(*s)*(g-h*(*tau));
}

void jacobi(double **a, int n, double *d, double **v)
{
  int j,iq,ip,i;
  int nrot=0;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  b=new double [n];
  z=new double [n];
  for (ip=0;ip<n;ip++) {
    v[ip][ip]=1.0;
    for (iq=ip+1;iq<n;iq++) v[ip][iq]=v[iq][ip]=0.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  for (i=0;i<50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) 
      for (iq=ip+1;iq<n;iq++)
	sm += fabs(a[ip][iq]);
    if (sm == 0.0) {
      delete [] z;
      delete [] b;
      break;
    }
    if (i < 3) tresh=0.2*sm/(n*n);
    else tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 3 && (fabs(d[ip])+g) == fabs(d[ip])
		  && (fabs(d[iq])+g) == fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((fabs(h)+g) == fabs(h)) t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
          for (j=0;j<ip;j++)
            ROTATE(a,j,ip,j,iq,&tau,&s);
	  for (j=ip+1;j<iq;j++)
	    ROTATE(a,ip,j,j,iq,&tau,&s);
	  for (j=iq+1;j<n;j++)
	    ROTATE(a,ip,j,iq,j,&tau,&s);
	  for (j=0;j<n;j++) 
	    ROTATE(v,j,ip,j,iq,&tau,&s);
	  ++(nrot);
	}
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  if(sm != 0.0) {
    printf("Too many iterations in routine jacobi %d ",nrot);
    if(DEBUG) {
      printf("original matrix\n");
      for(i=0;i<n;i++) {
	for(j=0;j<n;j++) printf("%8.4e ",a[i][j]);
	printf("\n");
      }
    }
    delete [] z;
    delete [] b;
  }
}

void eigsrt(double *d, double **v, int n)
{
  int k,j,i;

  for (i=0;i<n-1;i++) {
    for (j=i+1;j<n;j++) {
      if (d[j] > d[i]) {
	SWAP(&d[j],&d[i]);
	for(k=0;k<n;k++) SWAP(&v[k][i],&v[k][j]);
      }
    }
  }
}

/* The following subroutines calculates the integration of weighting    */
/* function for entropy from 0+ to 1/2 fmin. The algorithm is based on  */
/* mathematica result:                                                  */
/* w[u_] := (u/ (Exp[u] - 1) - Log[1 - Exp[-u]] )                       */
/* Integrate[  w[u], {u, 0, x}] =                                       */
/* (Pi^2)/2 -x^2 -x Log[1-Exp[-x]] + 2 PolyLog[2,Exp[-x]])              */

double scweighting(double upper)
{
  if(upper==0 || isnan(upper)) return 0;
  return (2*upper-upper*log(upper))/upper;
}

/* ---------------------------------------------------------------------- */

double sqweighting(double upper)
{
  if(upper==0 || isnan(upper)) return 0;
  double pi;
  double wsq;
  pi=3.1415926535897932385;

  wsq= pi*pi/3.0 -upper*upper +upper*log(-1+exp(upper))
     - 2*polylog(2.0,exp(-upper));

  return wsq/upper;
}

/* ---------------------------------------------------------------------- */

double polylog(double n, double z)
{
  int k;
  double sum,sum_old;

  k=1;
  sum=0.0;
  sum_old=0.0;

  do 
  {
    sum_old=sum;
    sum+= pow(z,k)/pow(k,n);
     k++;
  } while (fabs(sum-sum_old)!=0);

 return sum;
}

/* ---------------------------------------------------------------------- */

double HSDF(double *pwr,int lenth,double fmin,int nmol,double fract_f)
{
  double hsdf,tmpg,tmps;
  int j;
  hsdf=tmpg=tmps=0;
  for(j=0;j<lenth;j++) {
    twoPT(&tmpg,&tmps,pwr[0],pwr[j],fmin*j,nmol,fract_f);
    if(j==0 || j==lenth-1) { hsdf += tmpg*fmin*0.5;}
    else { hsdf += tmpg*fmin; }
  }
  return hsdf;
}

/* ---------------------------------------------------------------------- */

double TTDF(double *pwr,int lenth,double fmin)
{
  double ttdf;
  int j;
  ttdf=0;
  for(j=0;j<lenth;j++) {
    if(j==0 || j==lenth-1) { ttdf += pwr[j]*fmin*0.5;}
    else { ttdf += pwr[j]*fmin; }
  }
  return ttdf;
}

/* ---------------------------------------------------------------------- */

void HSweighting(double *wep,double *wap,double *wcvp,double *wsehs,double *wsp,double *wspd,double y,double mass,double nmol,double transT,double rotT,double volume,double *wsr,double *war,double *rT,double rs)
{
  *wep=*wap=*wcvp=*wsehs=*wsp=*wspd=*wsr=*war=0;
  //Ideal gas
  if(y<=0.74) { //packing fraction too large, ignore hard sphere correction
    *wsp=  5.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30/(nmol)); // indistinguishable particles 
    *wspd= 3.0/2.0 + log(pow(2.0*PI*mass*1e-3/Na*kb*transT/h/h,1.5)*volume*1e-30); // distinguishable particles 

    //Carnahan-Starling
    *wsehs=log((1.0+y+y*y-y*y*y)/pow(1.0-y,3.0))+y*(3.0*y-4.0)/pow(1.0-y,2.0);

    //Hard sphere gas
    *wsp  = (*wsp  + *wsehs)/3.0;
    *wspd = (*wspd + *wsehs)/3.0;
    *wep  = 0.5;
    *wap  = *wep-*wsp; /*wa=we-ws, no need of factor T in front of ws*/
    *wcvp = 0.5;

     //rigid rotor
    if(rT[0]>0.0001) { //for none-atomic molecules
      if(rT[2]<0) { //linear
        *wsr  = (1.0 + log(rotT/sqrt(rT[0]*rT[1])/rs))/2.0;
      } else { //non-linear
        *wsr  = (3.0/2.0 + log(sqrt(PI*rotT*rotT*rotT/(rT[0]*rT[1]*rT[2]))/rs))/3.0;
      }

      *war  = *wep-*wsr;
    } else *wsr = *war = 0;
  }

}

/* ---------------------------------------------------------------------- */

void twoPT(double *tmpg,double *tmps,double s0,double sv,double v,int nmol,double fract_f)
{
  if(fract_f==0) { *tmps=sv; *tmpg=0; }
  else if(v==0) { *tmpg=s0; *tmps=0; }
  else {
    *tmpg=s0/(1.0+pow(PI*s0*v/(6.0*nmol*fract_f),2.0));
    //*tmpg=3*fract_f*nmol*(s0/(1.0+pow(PI*s0*v/2,2.0)));
    if(*tmpg>sv) *tmpg=sv;
    *tmps=sv-*tmpg;
  }
}

/*determine fraction factor f from constant K*/

double search2PT(double K)
{
  double P,fold,fnew,dPdf,tol;
  int count;

  fold=0.0;
  fnew=0.7293*pow(K,0.5727); /*initial guess*/
  if(fnew>0.5) fnew=0.5;
  tol=1e-10;
  count=0;

  while( fabs(fnew-fold)>tol && count<999) {
    dPdf=0.0;
    fold=fnew;
    P= 2.0*pow(K,-4.5)*pow(fnew,7.5)-6.0*pow(K,-3.0)*pow(fnew,5.0)-pow(K,-1.5)*pow(fnew,3.5)
       + 6.0*pow(K,-1.5)*pow(fnew,2.5)+2.0*fnew-2;
    dPdf = 15.0*pow(K,-4.5)*pow(fnew,6.5)-30.0*pow(K,-3.0)*pow(fnew,4.0)-3.5*pow(K,-1.5)*pow(fnew,2.5)
           + 15.0*pow(K,-1.5)*pow(fnew,1.5)+2.0;
    fnew=fold-(P)/dPdf;
    count++;
  }
  return fnew;
}

/* ---------------------------------------------------------------------- */

void copyright(ostream *outf)
{
  char null[1024];

  strcpy(null,"**********************************************************************************");
  *outf<<null<<endl;
  strcpy(null,"*                     Two-Phase Thermodynamics (2PT) Program                     *");
  *outf<<null<<endl;
  strcpy(null,"*                          Shiang-Tai Lin (stlin@ntw.edu.tw)                     *");
  *outf<<null<<endl;
  strcpy(null,"* Department of Chemical Engineering, National Tiawan University, Taipei, Taiwan *");
  *outf<<null<<endl;
  strcpy(null,"*                   Prabal K Maiti (maiti@physics.iisc.ernet.in)                 *");
  *outf<<null<<endl;
  strcpy(null,"*   Department of Physics, Indian Institute of Science, Bangalore, India 560012  *");
  *outf<<null<<endl;
  strcpy(null,"*                          Tod A Pascal (tpascal@wag.caltech.edu)                *");
  *outf<<null<<endl;
  strcpy(null,"*                                     and                                        *");
  *outf<<null<<endl;
  strcpy(null,"*                        William A Goddard III (wag@wag.caltech.edu)             *");
  *outf<<null<<endl;
  strcpy(null,"*         Materials and Process Simulation Center, Caltech, Pasadena, CA USA     *");
  *outf<<null<<endl;
  strcpy(null,"*                                Copyright (c) 2010                              *");
  *outf<<null<<endl;
  strcpy(null,"*                                All rights Reserved                             *");
  *outf<<null<<endl;
  strcpy(null,"*                             *************************                          *");
  *outf<<null<<endl;
  strcpy(null,"*                                      Cite:                                     *");
  *outf<<null<<endl;
  strcpy(null,"* S.T. Lin, M. Blanco and W.A. Goddard, J. Chem. Phys., 2003, 119, 11792-11805   *");
  *outf<<null<<endl;
  strcpy(null,"* S.T. Lin, P.K. Maiti and W.A. Goddard, J. Phys. Chem. B., 2010, 114, 8191-8198 *");
  *outf<<null<<endl;
  strcpy(null,"* T.A. Pascal, S.T. Lin and W.A. Goddard, PCCP, 2011, 13 (1), 169 - 181          *");
  *outf<<null<<endl;
  strcpy(null,"**********************************************************************************");
  *outf<<null<<endl;
}
