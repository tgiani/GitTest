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

//one dimensional LO and NLO integrals. 
//no intrinsic component. bpdf is 0 below mb 

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
    double Asb = pow(as,2)*(L*L*asb(z,2,2) + L*asb(z,2,1)+asb(z,2,0));
    double fg  = (double) p_pdf->xfxQ(21,kin.x/z,kin.Q);
    double fsigma = fs(kin.x/z);
    
    return Agb*fg + Asb*fsigma;
  }
  // if Q<mb the bpdf is zero
  return 0.;
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

  bpdf = inte;
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
  APFEL::LHAPDFgrid(0,Qin,"NNPDF31_b_expanded__NLL");
 
  return 0;
}











