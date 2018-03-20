#include "iostream"
#include "iomanip"
#include "cmath"
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace APFEL;
using namespace LHAPDF;

LHAPDF::PDF * p_pdf;
double mb = 4.5;
double Q0 = 2;




  extern "C" void externalsetapfel_(const double& x, const double& Q, double *xf)
{
  
  for (int i=0; i<13; ++i){
    xf[i] = (double) p_pdf->xfxQ(i-6,x,Q);    // light flavour and charm PDFs are those of the input set, which for that flavours has NNPDF31
  }

  // for the bottom quark, if Q<mb we take NNPDF31_bintrinsic, if Q>mb we evolve NNPDF31_bintrinsic in 5FS
  double bpdf(0.);
  if(Q > mb){

       double Qa(mb), Qb(Q), eps = 1e-10;
       APFEL::SetPerturbativeOrder(2);
       APFEL::SetPDFSet("NNPDF31_bintrinsic_NLL.LHgrid");
       APFEL::InitializeAPFEL(); 
       APFEL::EvolveAPFEL(Q0,Qa);
       bpdf = APFEL::xPDF(5,x);
   
      }

      if(Q <= mb){

       bpdf = (double) p_pdf->xfxQ(5,x,Q);;
       
} 
  xf[1] = xf[11] = bpdf;
}



int main()
{

  APFEL::SetTheory("QCD");
  APFEL::SetPerturbativeOrder(2);
  APFEL::SetPDFSet("external");  
  APFEL::SetMaxFlavourAlpha(5);
  APFEL::SetMaxFlavourPDFs(5);
  APFEL::SetVFNS();                // apfel evolution will be used only when Q>mb for bpdf, therefore only with 5 active flavours

  p_pdf = LHAPDF::mkPDF("NNPDF31_nnlo_as_0118", 0);
  APFEL::InitializeAPFEL();

  APFEL::SetPoleMasses(p_pdf->quarkMass(4),
		       p_pdf->quarkMass(5),
		       p_pdf->quarkMass(6));

  APFEL::SetAlphaQCDRef(p_pdf->alphasQ(91.1987),91.1876);
  double Qin(-1.);
  //APFEL::InitializeAPFEL();
  APFEL::SetLHgridParameters(100,50,1e-9,1e-1,1,50,Qin*Qin,1e10);
  APFEL::LHAPDFgrid(0,Qin,"NNPDF31_test_NLL");
 
  return 0;

}
