#ifndef bpdf_H
#define bpdf_H
#include "LHAPDF/LHAPDF.h"

struct para {
  double Q0, mub, Q, x;
  LHAPDF::PDF * pdf;
};


double integrand_lo(double z, void *p);
double integrand_nlo(double z, void *p);
double intrinsic_component(double z[], size_t dim, void *p);
double integral_ic(const struct para p);
double integral_lo(const struct para p);
double integral_nlo(const struct para p);

#endif
