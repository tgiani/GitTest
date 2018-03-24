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

// compute the overall normalization A of fb^(4)(x) requiring that int_0^1  x*(fb^(4)(x) + fbbar^(4)(x)) = 1e-4


LHAPDF::PDF * p_pdf;

// definition of the integrand  x*fg0*Agb0
double integrand(double z[], size_t dim, void *p)
{
 
  double m2b   = p_pdf->alphaS().quarkMass(5)*p_pdf->alphaS().quarkMass(5),
         Q0    = 2,
         as0   = p_pdf->alphasQ(Q0)/(2.*M_PI),
         L0    = log(m2b/Q0*Q0);

  if(z[0]>z[1]){

    double Agb0 = as0*L0*agb(z[0],1,1);  // LO intrinsic component
    double fg0  = (double) p_pdf->xfxQ(21,z[1]/z[0],Q0);
    return z[1]*Agb0*fg0;

  }
  else{

    return 0;

  }

}


double integral()
{
  size_t calls = 200000;
  double inte(0.),w2(0.);
  size_t d=2;
  double xmin[d], xmax[d];
  for(size_t i(0); i<d; ++i){
    xmin[i] = 0.; xmax[i] =1.;
  }
  gsl_monte_function F;
  F.f = &integrand;
  F.dim = d;
  gsl_monte_vegas_state* state = gsl_monte_vegas_alloc(d);
  gsl_monte_vegas_init(state);
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus);
  gsl_monte_vegas_integrate(&F,xmin,xmax,d,calls,rng,state,&inte,&w2);
  gsl_monte_vegas_integrate(&F,xmin,xmax,d,calls,rng,state,&inte,&w2);
  gsl_rng_free(rng);
  gsl_monte_vegas_free(state);
  return 2*inte;
}


int main(){
 p_pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);
 double res = integral();
 cout << res << endl; 
 cout << "A :" <<  1e-4/res << endl;
 cout << "mb :" << p_pdf->alphaS().quarkMass(5) << endl;
 return 0;
}



