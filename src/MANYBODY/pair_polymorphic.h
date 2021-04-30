/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(polymorphic,PairPolymorphic)

#else

#ifndef LMP_PAIR_POLYMORPHIC_H
#define LMP_PAIR_POLYMORPHIC_H

#include "pair.h"
#include <cmath>

namespace LAMMPS_NS {

class PairPolymorphic : public Pair {

  public:

  PairPolymorphic(class LAMMPS *);
  virtual ~PairPolymorphic();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);

  protected:

  class tabularFunction {

    public:

    tabularFunction() {
      size = 0;
      xmin = 0.0;
      xmax = 0.0;
      xmaxsq = xmax*xmax;
      vmax = 0.0;
      xs = nullptr;
      ys = nullptr;
      ys1 = nullptr;
      ys2 = nullptr;
      ys3 = nullptr;
      ys4 = nullptr;
      ys5 = nullptr;
      ys6 = nullptr;
    }
    tabularFunction(int n) {
      size = n;
      xmin = 0.0;
      xmax = 0.0;
      xmaxsq = xmax*xmax;
      xs = new double[n];
      ys = new double[n];
      ys1 = new double[n];
      ys2 = new double[n];
      ys3 = new double[n];
      ys4 = new double[n];
      ys5 = new double[n];
      ys6 = new double[n];
    }
    tabularFunction(int n, double x1, double x2) {
      size = n;
      xmin = x1;
      xmax = x2;
      xmaxsq = xmax*xmax;
      xs = new double[n];
      ys = new double[n];
      ys1 = new double[n];
      ys2 = new double[n];
      ys3 = new double[n];
      ys4 = new double[n];
      ys5 = new double[n];
      ys6 = new double[n];
    }
    virtual ~tabularFunction() {
      delete [] xs;
      delete [] ys;
      delete [] ys1;
      delete [] ys2;
      delete [] ys3;
      delete [] ys4;
      delete [] ys5;
      delete [] ys6;
    }
    void set_xrange(double x1, double x2) {
      xmin = x1;
      xmax = x2;
      xmaxsq = xmax*xmax;
    }
    void set_values(int n, double x1, double x2, double * values, double epsilon)
    {
      int shrink = 1;
      int ilo,ihi;
      double vlo,vhi;
      ilo = 0;
      ihi = n-1;
      for (int i = 0; i < n; i++) {
        if (fabs(values[i]) <= epsilon) {
          ilo = i;
        } else {
          break;
        }
      }
      for (int i = n-1; i >= 0; i--) {
        if (fabs(values[i]) <= epsilon) {
          ihi = i;
        } else {
          break;
        }
      }
      if (ihi < ilo) ihi = ilo;
      vlo = values[ilo];
      vhi = values[ilo];
      for (int i = ilo; i <= ihi; i++) {
        if (vlo > values[i]) vlo = values[i];
        if (vhi < values[i]) vhi = values[i];
      }

//    shrink (remove near zero points) reduces cutoff radius, and therefore computational cost
//    do not shrink when x2 < 1.1 (angular function) or x2 > 20.0 (non-radial function)
      if (x2 < 1.1 || x2 > 20.0) {
        shrink = 0;
      }
//    do not shrink when when list is abnormally small
      if (ihi - ilo < 50) {
        shrink = 0;
      }
//    shrink if it is a constant
      if (vhi - vlo <= epsilon) {
//        shrink = 1;
      }

      if (shrink == 0) {
        ilo = 0;
        ihi = n-1;
      }
      xmin = x1 + (x2-x1)/(n -1)*ilo;
      xmax = xmin + (x2-x1)/(n -1)*(ihi-ilo);
      xmaxsq = xmax*xmax;
      n = ihi - ilo + 1;
      resize(n);
      for (int i = ilo; i <= ihi; i++) {
        ys[i-ilo] = values[i];
      }
      initialize();
    }
    void value(double x, double &y, int ny, double &y1, int ny1)
    {
      double ps = (x - xmin) * rdx;
      int ks = ps + 0.5;
      if (ks > size-1) ks = size-1;
      if (ks < 0 ) ks = 0;
      ps = ps - ks;
      if (ny) y = ((ys3[ks]*ps + ys2[ks])*ps + ys1[ks])*ps + ys[ks];
      if (ny1) y1 = (ys6[ks]*ps + ys5[ks])*ps + ys4[ks];
    }
    void print_value()
    {
      printf("%d %f %f %f \n",size,xmin,xmax,rdx);
      printf(" \n");
      for (int i = 0; i < size; i++) {
        printf("%f %f \n",xs[i],ys[i]);
      }
    }
    double get_xmin() {
      return xmin;
    }
    double get_xmax() {
      return xmax;
    }
    double get_xmaxsq() {
      return xmaxsq;
    }
    double get_rdx() {
      return rdx;
    }
    double get_vmax() {
      return vmax;
    }

    protected:

    void resize(int n) {
      if (n != size) {
        size = n;
        delete [] xs;
        xs = new double[n];
        delete [] ys;
        ys = new double[n];
        delete [] ys1;
        ys1 = new double[n];
        delete [] ys2;
        ys2 = new double[n];
        delete [] ys3;
        ys3 = new double[n];
        delete [] ys4;
        ys4 = new double[n];
        delete [] ys5;
        ys5 = new double[n];
        delete [] ys6;
        ys6 = new double[n];
      }
    }
    void initialize() {
      int n = size;
      rdx = (xmax-xmin)/(n-1.0);
      vmax = 0.0;
      for (int i = 0; i < n; i++) {
        if (fabs(ys[i]) > vmax) vmax = fabs(ys[i]);
      }
      for (int i = 0; i < n; i++) {
        xs[i] = xmin+i*rdx;
      }
      rdx = 1.0 / rdx;
      ys1[0] = ys[1] - ys[0];
      ys1[1] = 0.5 * (ys[2] - ys[0]);
      ys1[n-2] = 0.5 * (ys[n-1] - ys[n-3]);
      ys1[n-1] = ys[n-1] - ys[n-2];
      for (int i = 2; i < n-2; i++) {
        ys1[i]=((ys[i-2]-ys[i+2])+ 8.0*(ys[i+1]-ys[i-1]))/12.0;
      }
      for (int i = 0; i < n-1; i++) {
        ys2[i]=3.0*(ys[i+1]-ys[i])-2.0*ys1[i]-ys1[i+1];
        ys3[i]=ys1[i]+ys1[i+1]-2.0*(ys[i+1]-ys[i]);
      }
      ys2[n-1]=0.0;
      ys3[n-1]=0.0;
      for (int i = 0; i < n; i++) {
        ys4[i]=ys1[i]*rdx;
        ys5[i]=2.0*ys2[i]*rdx;
        ys6[i]=3.0*ys3[i]*rdx;
      }
    }
    int size;
    double xmin,xmax,xmaxsq,rdx,vmax;
    double *ys, *ys1, *ys2, *ys3, *ys4, *ys5, *ys6;
    double *xs;
  };

  struct PairParameters {
    double cut;
    double cutsq;
    double xi;
    class tabularFunction *U;
    class tabularFunction *V;
    class tabularFunction *W;
    class tabularFunction *F;
    PairParameters() {
      cut = 0.0;
      cutsq = 0.0;
      xi =  1.0;
      U = nullptr;
      V = nullptr;
      W = nullptr;
      F = nullptr;
    };
    ~PairParameters() {
      delete U;
      delete V;
      delete W;
      delete F;
    }
  };
  struct TripletParameters {
    class tabularFunction *P;
    class tabularFunction *G;
    TripletParameters() {
      P = nullptr;
      G = nullptr;
    };
    ~TripletParameters() {
      delete P;
      delete G;
    }
  };

  double epsilon;
  int eta;
  int nx,nr,ng; // table sizes
  double maxX;

  // parameter sets
  PairParameters    * pairParameters;    // for I-J interaction
  TripletParameters * tripletParameters; // for I-J-K interaction

  int neighsize,numneighV,numneighW,numneighW1;
  int *firstneighV,*firstneighW,*firstneighW1;
  double *delxV,*delyV,*delzV,*drV;
  double *delxW,*delyW,*delzW,*drW;

  double cutmax;                // max cutoff for all elements
  double cutmaxsq;
  int npair,ntriple;
  int *match;

  void allocate();

  virtual void read_file(char *);
  void setup_params();
  void write_tables(int);

  void attractive(PairParameters *, PairParameters *, TripletParameters *,
                  double, double, double, double *, double *, double *,
                  double *, double *);

  void ters_zetaterm_d(double, double *, double, double *, double, double *,
                       double *, double *, PairParameters *, PairParameters *,
                       TripletParameters *);
  void costheta_d(double *, double, double *, double,
                  double *, double *, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style polymorphic requires atom IDs

This is a requirement to use the polymorphic potential.

E: Pair style polymorphic requires newton pair on

See the newton command.  This is a restriction to use the polymorphic
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open polymorphic potential file %s

The specified polymorphic potential file cannot be opened.  Check that
the path and name are correct.

E: Incorrect number of elements in potential file

Self-explanatory.

E: Element not defined in potential file

The specified element is not in the potential file.

E: Potential file incompatible with this pair style version

UNDOCUMENTED

E: Error reading potential file header

UNDOCUMENTED

*/
