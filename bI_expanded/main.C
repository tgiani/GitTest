#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <getopt.h>
#include "helpers.h"

using namespace LHAPDF;
using namespace std;


LHAPDF::PDF * p_pdf;

struct Kinetics {
  double x;
  double Q;
};

Kinetics kin;

double fs(double z)
{
  double res = 0;
  for(size_t i(1); i<5 ; ++i){
    res+=(double) (p_pdf->xfxQ(i,z,kin.Q)+p_pdf->xfxQ(-i,z,kin.Q));
  }
  
  return res;
}

//one dimensional LO and NLO integrals

double integrand(double z, void * p)
{
  double Q = kin.Q,
         as = p_pdf->alphasQ(Q)/(2.*M_PI),
         m2b   = p_pdf->alphaS().quarkMass(5)*p_pdf->alphaS().quarkMass(5),
         Q0 = 2,
         L     = log(Q*Q/m2b),
         as0   = p_pdf->alphasQ(Q0)/(2.*M_PI),
         L0    = log(m2b/Q0*Q0);

  if(z>kin.x && m2b<Q*Q){

    double Agb = as*L*agb(z,1,1) 
      + pow(as,2)*(L*L*agb(z,2,2) + L*agb(z,2,1)+agb(z,2,0));
    double Agb0 = as0*L0*agb(z,1,1);  // LO intrinsic component
    double Asb = pow(as,2)*(L*L*asb(z,2,2) + L*asb(z,2,1)+asb(z,2,0));
    double Abb_ep = - 2.*abb_ep(kin.x/z,L)*Pqg(1.)*p_pdf->xfxQ(21,z,Q0)/z *as*as0*L0;  // NLO intrinsic componenet
    double fg  = (double) p_pdf->xfxQ(21,kin.x/z,kin.Q);
    double fg0  = (double) p_pdf->xfxQ(21,kin.x/z,Q0);
    double fsigma = fs(kin.x/z);
    
    return Agb*fg + Asb*fsigma + Agb0*fg0 + Abb_ep;
  }
  // if Q<mb I keep the constant part due to thre intrinsic componenet
  if(z>kin.x && m2b>=Q*Q){

    double Agb0 = as0*L0*agb(z,1,1);  // LO intrinsic component
    double Abb_ep = - 2.*abb_ep(kin.x/z,L)*Pqg(1.)*p_pdf->xfxQ(21,z,Q0)/z *as*as0*L0;  // NLO intrinsic componenet
    double fg0  = (double) p_pdf->xfxQ(21,kin.x/z,Q0);

    return Agb0*fg0 + Abb_ep;
  }

  return 0.;
}

//two dimensional integral coming from the intrinsic componenet

double intrinsic_component(double z[], size_t dim, void *p)
{

 double Q = kin.Q;
 double as = p_pdf->alphasQ(Q)/(2.*M_PI),
       m2b   = p_pdf->alphaS().quarkMass(5)*p_pdf->alphaS().quarkMass(5),
       Q0 = 2,
       L     = log(Q*Q/m2b),
       as0   = p_pdf->alphasQ(Q0)/(2.*M_PI),
       L0    = log(m2b/Q0*Q0);

  double res(0.);
  double tau, jac(1.), xi, ximax;
  double zz,tt;

  tau = kin.x*pow(1./kin.x,z[0]);
  jac *= log(1./kin.x)*tau;
  
  ximax = -1./2.*log(tau);
  xi    = -ximax + 2.*ximax*z[1];
  jac  *= 2.*ximax;

  tt = sqrt(tau)*exp(xi);
  zz = sqrt(tau)*exp(-xi);
  if(Q*Q>m2b){

    if(1.-zz > 1.e-5){
      res = jac*abb(zz,L)*(Pqg(kin.x/zz/tt) - zz*Pqg(kin.x/tt))/tau;  //#### pqg-->Pqg
    }
    return 2.*res*p_pdf->xfxQ(21,tt,Q0)*as*as0*L0;                    //#### t-->tt
  }
  // if Q<mb I keep the constant part due to thre intrinsic componenet
  else{

    if(1.-zz > 1.e-5){
      res = jac*abb(zz,0)*(Pqg(kin.x/zz/tt) - zz*Pqg(kin.x/tt))/tau;  //#### pqg-->Pqg
    }
    return 2.*res*p_pdf->xfxQ(21,tt,Q0)*as0*as0*L0;                   //#### t-->tt 
  }
    
}                


extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  kin.x = x;
  kin.Q = Q;

  
  for (int i=0; i<13; ++i){
    xf[i] = (double) p_pdf->xfxQ(i-6,x,Q);
  }

  double bpdf(0.);;
  double inte(0.),res(0.),err(0);
  size_t neval(0);
  gsl_integration_cquad_workspace * state = gsl_integration_cquad_workspace_alloc(100);
  gsl_function F;
  F.function = &integrand;
  gsl_integration_cquad(&F,kin.x,p_pdf->xMax(),1.e-12,1.e-12,state,&inte,&err,&neval);
  gsl_integration_cquad_workspace_free(state);
  res = inte;

  size_t calls = 200000;
  double inte1(0.),w2(0.);
  size_t d=2;
 
  double xmin[d], xmax[d];
  for(size_t i(0); i<d; ++i){
    xmin[i] = 0.; xmax[i] =1.;                         //#### why now the integrals are performed both between 0 and 1 ?
  }
 
  gsl_monte_function G;
  G.f = &intrinsic_component;
  G.dim = d;
  G.params = &kin;
  gsl_monte_vegas_state* state1 = gsl_monte_vegas_alloc(d);
  gsl_monte_vegas_init(state1);
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&G,xmin,xmax,d,calls,rng,state1,&inte1,&w2);
  gsl_monte_vegas_integrate(&G,xmin,xmax,d,calls,rng,state1,&inte1,&w2);
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state1);
  
  bpdf = inte + inte1 ;
  xf[1] = xf[11] = bpdf;
}


int main()
{
  APFEL::SetTheory("QCD");
  APFEL::SetPerturbativeOrder(2);
  APFEL::SetPDFSet("external");
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);

  p_pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);
  APFEL::SetPoleMasses(p_pdf->quarkMass(4),
		       p_pdf->quarkMass(5),
		       p_pdf->quarkMass(6));

  APFEL::SetAlphaQCDRef(p_pdf->alphasQ(91.1987),91.1876);
  double Qin(-1.);
  APFEL::SetVFNS();
  //APFEL::InitializeAPFEL();
  //APFEL::SetLHgridParameters(100,50,1e-9,1e-1,1,50,Qin*Qin,1e10);
  APFEL::LHAPDFgrid(0,Qin,"NNPDF31_bI_expanded__NLL");
 
  return 0;
}











