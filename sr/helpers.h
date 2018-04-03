#ifndef HELPERS_H
#define HELPERS_H
#include "LHAPDF/LHAPDF.h"
extern "C" {
  double wgplg_(int *N, int*P, double *x);
}
double dilog(double z);
double trilog(double z);
double Pqg(const double x);
double agb(const double z, const int fo, const int lo);
double asb(const double z, const int fo, const int lo);
double pqq(const double z);
double abb(const double z, const double L);
double abb_ep(const double x, const double L);
double sigma(const double x, const double Q, const LHAPDF::PDF * pdf);
#endif
