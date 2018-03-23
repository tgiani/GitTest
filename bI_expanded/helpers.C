#include "helpers.h"
#include "LHAPDF/LHAPDF.h"

#define ZETA2 1.6449340668482264365
#define ZETA3 1.2020569031595942854



/* di2->dilog, di3->trilog */


double dilog(double z) // PolyLog[2,z] or HPL[{2},z]==HPL[{0,1},z]
{

  const double pi2o6 = 1.644934066848226; // zeta2

  static double cof_odd[9] =
    { 1.,0.027777777777777776,-0.0002777777777777778,4.72411186696901e-6,
      -9.185773074661964e-8,1.8978869988971e-9,4.0647616451442256e-11,
      8.921691020456452e-13,-1.9939295860721074e-14 };

  if ( z==0. ) return 0.;
  else if ( z==1. ) return pi2o6;
  else if ( z==-1. ) return -pi2o6/2.;
  else if ( -1<z && z<=0.5 ) {
    double u = -log(1.-z);
    double t = u*u, f = 0.;
    int i; for (i=8; i>=0; i--) f = f*t+cof_odd[i];
    return f*u - t/4.;
  }

  if ( z<-1. ) {
    double t = 1./z;
    double lnnt = log(-t);
    //return -pi2o6 - li2(t) - lnnt*lnnt/2.;
    return -pi2o6 - dilog(t) - lnnt*lnnt/2.;
  }
  else if ( 0.5<z && z<1. ) {
    double t = 1.-z;
   // return pi2o6 - log(t)*log(1-t) - li2(t);
   return pi2o6 - log(t)*log(1-t) - dilog(t);
  }
  else if ( 1.<z && z<2. ) {
    double t = 1.-1./z;
    double lnmt = log(1.-t), lnt = log(t);
    // omit imaginary part of pi*log(z)*sign(Im(z))
    // return pi2o6 - lnmt*lnmt/2. + lnmt*lnt + li2(t);
    return pi2o6 - lnmt*lnmt/2. + lnmt*lnt + dilog(t); 
  }
  else if ( z>=2. ) {
    double t = 1./z;
    double lnt = log(t);
    // omit imaginary part of pi*log(z)*sign(Im(z))
    //return 2.*pi2o6 - li2(t) - lnt*lnt/2.;
    return 2.*pi2o6 - dilog(t) - lnt*lnt/2.; 
  }

}

double trilog(double z) // PolyLog[3,z] or HPL[{3},z]==HPL[{0,0,0,1},z]
{

  const double pi2o6 = 1.644934066848226; // zeta2
  const double zeta3 = 1.202056903159594;

  static double cof[17] =
    { 1.,-0.375,0.07870370370370369,-0.008680555555555556,0.00012962962962962963,
      0.00008101851851851852,-3.4193571608537595e-6,-1.328656462585034e-6,
      8.660871756109852e-8,2.52608759553204e-8,-2.1446944683640645e-9,
      -5.140110622012978e-10,5.249582114600829e-11,1.0887754406636318e-11,
      -1.2779396094493697e-12,-2.369824177308745e-13,3.104357887965462e-14 };
  static double cofm_even[9] = 
    { 1.202056903159594,0.75,-0.003472222222222222,0.000011574074074074073,
      -9.841899722852104e-8,1.1482216343327454e-9,-1.5815724990809165e-11,
      2.4195009792525154e-13,-3.982897776989488e-15 };

  if ( z==0. ) return 0.;
  else if ( z==1. ) return zeta3;
  else if ( z==-1. ) return -3.*zeta3/4.;
  else if ( -1.<z && z<0.5 ) {
    double u = -log(1.-z);
    double f = 0.;
    int i; for (i=16; i>=0; i--) f = f*u+cof[i];
    return f*u;
  }
  else if ( 0.5<=z && z<1. ) {
    double u = -log(z);
    double t = u*u, f = 0.;
    int i; for (i=8; i>=0; i--) f = f*t+cofm_even[i];
    return f - pi2o6*u + u*t/12. - t*log(u)/2.;
  }

  if ( z<-1. ) {
    double t = 1./z;
    double lnnt = log(-t);
   // return li3(t) + pi2o6*lnnt + lnnt*lnnt*lnnt/6.;
   return trilog(t) + pi2o6*lnnt + lnnt*lnnt*lnnt/6.;
  }
  else if ( z>1. ) {
    double t = 1./z;
    double lnt = log(t);
    // omit imaginary part of pi/2*log(z)^2*sign(Im(z))
    //return -2.*pi2o6*lnt + lnt*lnt*lnt/6. + li3(t);
    return -2.*pi2o6*lnt + lnt*lnt*lnt/6. + trilog(t); 
  }

}

//Matching coefficeint


double Pqg(const double x)                       //### pqg --> Pqg
{
  return 0.5*(x*x+pow((1.0-x),2));
}


double agb(const double z, const int fo, const int lo)
{
  double CF(4./3.), CA(3.), TR(0.5), lz(log(z)),lIz(log(1.-z)),
      z2(z*z);
  double res(0.), cftr(0.), catr(0.), tr2(0.);
  if(fo==1 && lo==1) res=Pqg(z);
  if(fo==2 && lo==0) {
      double Iz(1. - z);
      int n(1),p(2);
      cftr = 2.*Pqg(z)*(8.*ZETA3 + 4./3.*pow(lIz,3) - 8.*lIz*dilog(Iz)
                        +8.*ZETA2*lz - 4.*lz*pow(lIz,2) + 2./3.*pow(lz,3)
                        - 8.*lz*dilog(1.-z) + 8.*trilog(1.-z) - 24.*wgplg_(&n,&p,&Iz))
          + z2*(-16.*ZETA2*lz + 4./3.*pow(lz,3) + 16.*lz*dilog(Iz)+32.*wgplg_(&n,&p,&Iz))
          - (4. + 96.*z - 64.*z2)*dilog(Iz) - (4.-48.*z + 40.*z2)*ZETA2 - (8.+48.*z - 24.*z2)*lz*lIz
          + (4.+8.*z- 12.*z2)*pow(lIz,2) - (1.+12.*z-20.*z2)*pow(lz,2)
          - (52.*z - 48.*z2)*lIz - (16.+18.*z+48.*z2)*lz
          + 26. - 82.*z + 80.*z2;

      double mz = -z;
      catr = 2.*Pqg(z)*(-4./3.*pow(lIz,3) + 8.*lIz*dilog(Iz) - 8.*trilog(Iz))
          +(1. + 2.*z + 2.*z2)*(-8.*ZETA2*log(1.+z) - 16.*log(1.+z)*dilog(-z)
                                -8.*lz*pow(log(1.+z),2) + 4.*pow(lz,2)*log(1.+z) + 8.*lz*dilog(-z)
                                -8.*trilog(-z) - 16.*wgplg_(&n,&p,&mz))
          + (16. + 64.*z)*(2.*wgplg_(&n,&p,&Iz)+lz*dilog(Iz))
          -(4./3. + 8./3.*z)*pow(lz,3)
          + (8.-32.*z + 16.*z2)*ZETA3
          -(16.+64.*z)*ZETA2*lz
          + (16.*z + 16.*z2)*(dilog(-z) +
                              lz*log(1.+z))
          + (32./3./z + 12. + 64.*z - 272./3.*z2)*dilog(Iz)
          - (12. + 48.*z - 260./3.*z2 + 32./3./z)*ZETA2
          - 4.*z2*lz*lIz
          - (2.+8.*z- 10.*z2)*pow(lIz,2)
          + (2. + 8.*z + 46./3.*z2)*pow(lz,2)
          + (4.+16.*z - 16.*z2)*lIz
          - (56./3. + 172./3.*z + 1600./9.*z2)*lz
          - 448./27./z - 4./3. - 628./3.*z + 6352./27.*z2;

      res = (CF*TR*cftr + CA*TR*catr)/8.;  // 8 is = 4 wrt Maria & 2 wrt to Buza...

    }
  if(fo==2 && lo==1) {
      cftr = 8.*Pqg(z)*(2.*lz*lIz - lIz*lIz + 2.*ZETA2)
          -(2. - 4.*z + 8.*z2)*lz*lz - 16.*z*(1.-z)*lIz
          -(6. - 8.*z + 16.*z2)*lz - 28. + 58.*z - 40.*z2;

      catr = (8. + 16.*z + 16.*z2)*(dilog(-z)+log(z)*log(1+z))
          + 8.*Pqg(z)*lIz*lIz
          + (4. + 8.*z)*lz*lz + 16.*z*ZETA2 + 16.*z*(1.-z)*lIz
          - (4. + 32.*z + 176.*z2/3.)*lz
          - 80./(9.*z) + 8. - 100.*z + 872.*z2/9.;

      res = - (CF*TR*cftr + CA*TR*catr)/4.;
    }
  if(fo==2 && lo==2){
      cftr = 8.*Pqg(z)*lIz - (2. - 4.*z + 8.*z2)*lz - (1. - 4.*z);

      catr = 8.*Pqg(z)*lIz + (4. + 16.*z)*lz
          + 8./(3.*z) + 2. + 16.*z - 62.*z2/3.;

      tr2 = -8.*(z2+(1.-z)*(1.-z))/3.;
      double tr(0.);
      tr = 8.*Pqg(z)/3.; // 8./3.*Pqg(z);

      res = (CF*TR*cftr - CA*TR*catr + TR*TR*tr2 + TR*tr)/4.;
    }
  return res;
}


double asb(const double z, const int fo, const int lo)
{
  double CF(4./3.), TR(0.5), lz(log(z)), z2(z*z);
  double res(0.), cftr(0.), lIz(log(1.-z));
  if(fo==2 && lo==0) {
      double a1(0.), a2(0.), a3(0.), a4(0.), a5(0.),
          a6(0.);
      double Iz(1. - z);
      int n(1),p(2);

      a1 = 32.*wgplg_(&n,&p,&Iz) + 16.*lz*dilog(Iz) - 16.*ZETA2*lz
          - 4./3.*pow(lz,3);
      a2 = 32./3./z + 8. -8.*z - 32./3.*z2;
      a3 = -a2;
      a4 = 2. + 10.*z+ 16./3.*z2;
      a5 = -(56./3.+88./3.*z + 448./9.*z2);
      a6 = -448./27./z - 4./3. - 124./3.*z + 1600./27.*z2;

      res =CF*TR*((1.+z)*a1 + dilog(Iz)*a2 + ZETA2*a3 + pow(lz,2)*a4
                  + lz*a5 + a6);
      res /= 8.;
//      res = 0.;

    }
  if(fo==2 && lo==2) {
      cftr = - 4.*(1. + z)*lz - 8./(3.*z) - 2. + 2.*z + 8.*z2/3. ;
      res = (CF*TR*cftr)/4.;
    }
  if(fo==2 && lo==1){
      cftr =- 4.*(1. + z)*lz*lz + (4. + 20.*z + 32.*z2/3.)*lz
          + 80./(9.*z) - 8. + 24.*z - 224.*z2/9.;

      res = (CF*TR*cftr)/4.;
    }
  return res;
}



//additional matching coefficient abb (series is in power of alpha/2Pi)

double pqq(const double z)
{
  if(1.-z < 1.e-5) return 0.; // to avoid instabilities with the plus dist
  return (1+sqrt(z))/(1.-z);
}

double abb(const double z, const double L)
{
  double res(0.),  CF(4./3.);
  res = CF*pqq(z)*(L-2.*log(1.-z)-1.); //Eq.(7) !!!
  return res;
}

double abb_ep(const double x, const double L) // end point contribution
{
  double  CF(4./3.);
  return CF*(-(x*(4 + L*(2 + x)))/2.
	     + log(1 - x)*(-1 - 2*L + 2*x + pow(x,2))
	     + 2*pow(log(1 - x),2));
}


/*
double sigma(const double x, const double Q, const LHAPDF::PDF * pdf)
{
  double CF(4./3.), res(0.);
  for(size_t i(1); i<5; ++i){
    res+= (pdf->xfxQ(i,x,Q)+pdf->xfxQ(-i,x,Q));
  }
  return res;
}
*/
