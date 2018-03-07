#include "iostream"
#include "iomanip"
#include "cmath"
#include "APFEL/APFEL.h"
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace APFEL;
using namespace LHAPDF;


double Q0=2;
double mb=4.5;

double xpdf(double x, double Q, int flavour){

  if(flavour==5 || flavour==-5)
          
      if(Q > mb){

       double Qa(mb), Qb(Q), eps = 1e-10;    
       APFEL::SetPerturbativeOrder(2);
       APFEL::SetFFNS(5);
       APFEL::SetPDFSet("NNPDF31_bintrinsic_NLL.LHgrid");
       APFEL::InitializeAPFEL();
       APFEL::EvolveAPFEL(Q0,Qa);
       return APFEL::xPDF(flavour,x);
   
      }

      if(Q <= mb){

       LHAPDF::PDF * below_mb_pdf = LHAPDF::mkPDF("NNPDF31_bintrinsic_NLL", 0); //pdf to use for Q<m
       return(below_mb_pdf->xfxQ(flavour,x,Q));
       
      } 


  else{

      if(Q > Q0){ 

       double Qa(mb), Qb(Q), eps = 1e-10;    
       APFEL::SetPerturbativeOrder(2);
       APFEL::SetVFNS;
       APFEL::SetPDFSet("NNPDF31_bintrinsic_NLL.LHgrid");
       APFEL::InitializeAPFEL();
       APFEL::EvolveAPFEL(Q0,Qa);
       return APFEL::xPDF(flavour,x);
   
      }

      if(Q <= Q0){

       LHAPDF::PDF * below_mb_pdf = LHAPDF::mkPDF("NNPDF31_bintrinsic_NLL", 0); 
       return(below_mb_pdf->xfxQ(flavour,x,Q));
       
      } 
            
  }

}


  
  

int main(){

 cout << "culo" << endl;
 return 0;
}
