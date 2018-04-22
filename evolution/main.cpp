#include <iostream>
#include <sstream>
#include <string.h>
#include <vector>
#include <list>
#include <fstream>
#include <algorithm>
#include <glob.h>
#include <iomanip>
#include <locale>
#include "stdio.h"
#include "stdlib.h"
#include <cassert>

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL libs
#include "APFEL/APFEL.h"

//gsl integration
#include <gsl/gsl_integration.h>

#include <unistd.h>

#define GetCurrentDir getcwd

string GetCurrentWorkingDir()
{
  char buff[FILENAME_MAX];
  GetCurrentDir( buff, FILENAME_MAX );
  string current_working_dir(buff);
  return current_working_dir;
}


const string set     = "NNPDF31_bI_expandedNORM_NLL";
LHAPDF::PDF * p_pdf  = LHAPDF::mkPDF(set,0);

const string st_set     = "NNPDF31_nnlo_as_0118";
LHAPDF::PDF * p_pdf_st  = LHAPDF::mkPDF(set,0);

const auto    pto    = p_pdf->orderQCD();
const auto    Qref   = p_pdf->info().get_entry_as<double>("MZ");
const auto    asref  = p_pdf->alphasQ(Qref);
const auto    muc    = p_pdf->quarkThreshold(4);
const auto    mc     = p_pdf->quarkMass(4);
const auto    mub    = p_pdf->quarkThreshold(5);
const auto    mb     = p_pdf->quarkMass(5);
const auto    mut    = p_pdf->quarkThreshold(6);
const auto    mt     = p_pdf->quarkMass(6);
const auto    Qmin   = p_pdf->qMin();
const auto    Qmax   = p_pdf->qMax();
const auto    xmin   = p_pdf->xMin();
const int     nrep   = 1;

const double threshold = 2.*mub;

struct inputs{ size_t i; double mu; };

struct inputs par;

double pdf_integrator(double x, void *p)
{
  return p_pdf->xfxQ(par.i,x,par.mu);
}

double st_integrator(double x, void *p)
{
  return p_pdf_st->xfxQ(par.i,x,par.mu);
}

double w(const size_t i, const double mu, const size_t mode = 0) 
{
  // returns int_0^1 dx x f_i (x,mu^2)
  par.i = i; par.mu = mu;
  double inte(0.),res(0.),err(0);
  size_t neval(0);
  gsl_integration_cquad_workspace * state = gsl_integration_cquad_workspace_alloc(100);
  gsl_function F;
  F.function = (mode)?&pdf_integrator:&st_integrator;
  gsl_integration_cquad(&F,1.e-12,p_pdf->xMax(),1.e-12,1.e-12,state,&inte,&err,&neval);
  gsl_integration_cquad_workspace_free(state);
  res = inte;
  return res;
}
double kfacx;
void init_kfacx(const double Q)
{
  kfacx = (w(5,Q,1)-w(5,Q)+w(-5,Q,1)-w(-5,Q) )/w(0,Q);
}

extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  for (int i=0; i<13; ++i)
    xf[i] = (double) p_pdf->xfxQ(i-6,x,Q);
  if(Q == threshold){
    //shift the gluon to have sum rules
    //not sure about the anti-b!!!
    xf[6] *= (1.-kfacx);
  }
}

int main(int argc, const char * argv[]) {
  std::cout << " mc : " << mc << ", muc : " << muc << "\n"
	    << " mb : " << mb << ", mub : " << mut << "\n"
	    << " mt : " << mt << ", mut : " << mut << "\n";

  init_kfacx(threshold);

  const string out_name      = "NNPDF31_bintrinsic_sum_rules";
  APFEL::EnableWelcomeMessage(false);
  APFEL::SetPDFSet("external");
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetTheory("QCD");
  APFEL::SetAlphaEvolution("expanded");
  APFEL::SetPDFEvolution("truncated");
  APFEL::SetPerturbativeOrder(pto);
  APFEL::SetAlphaQCDRef(asref, Qref);
  APFEL::SetVFNS();
  APFEL::SetPoleMasses(mc, mb, mt);
  APFEL::SetMassMatchingScales(1,1,1);
  APFEL::SetLHgridParameters(150, 75, xmin, 0.1, 1, 50, Qmin*Qmin, Qmax*Qmax);
  APFEL::LockGrids(true);
  APFEL::LHAPDFgrid(nrep-1, threshold, out_name);
  APFEL::CleanUp();
  
  return 0;
}
