#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "GeneralGrid.h"
#include <iostream>
#include <cstdlib>


/// The CubicSpline class is a third-order spline representation of a
/// function.  It stores a pointer to a grid and the values of the
/// function and its second derivative at the points defined by the
/// grid.
class CubicSpline {

 private:
  /// This flag records whether or not the stored second derivatives
  /// are in sync with the function values.  It is used to determine
  /// whether the second derivatives need recomputation.
  bool UpToDate;
  /// The function values on the grid points.
  vector<double> y;
  /// The second derivatives of the function
  vector<double> d2y;
  /// The values of the derivative of the represented function on the
  /// boundary.  If each value is greater that 1e30, we compute
  /// bondary conditions assuming that the second derivative is zero at
  /// that boundary.
  double StartDeriv, EndDeriv;

 public:
  GeneralGrid grid;
  /// Returns the interpolated value.
  inline int size() { return grid.NumPoints(); }
  inline double operator()(double x);
  /// Returns the interpolated first derivative.
  inline double Deriv(double x);
  /// Returns the interpolated second derivative.
  inline double Deriv2(double x);
  /// Returns the interpolated third derivative.
  inline double Deriv3(double x);
  /// Recompute the second derivatives from the function values
  void Update();
  /// Initialize the cubic spline.  See notes about start and end
  /// deriv above.
  inline void Init(GeneralGrid &newGrid, vector<double> yvals,
   double startderiv, double endderiv) {
   StartDeriv = startderiv;
   EndDeriv   = endderiv;
   if (newGrid.NumPoints() != yvals.size()) {
    cerr << "Size mismatch in CubicSpline.\n";
    cerr << "Grid Points = " << newGrid.NumPoints() << endl;
    cerr << "Y points    = " << yvals.size() << endl;
    abort();
   }
   grid = newGrid;
   y.resize(grid.NumPoints());
   d2y.resize(grid.NumPoints());
   y = yvals;
   Update();
  }
  /// Simplified form which assumes that the second derivative at both
  /// boundaries are zero.
  inline void Init (GeneralGrid &newGrid, vector<double> &yvals) {
   Init (newGrid, yvals, 5.0e30, 5.0e30);
  }
  /// Simplified constructor.
  inline CubicSpline (GeneralGrid &newGrid, vector<double> &yvals) {
   StartDeriv = EndDeriv = 5.0e30;
   Init (newGrid, yvals, 5.0e30, 5.0e30);
  }
  /// Full constructor.
  inline CubicSpline (GeneralGrid &newGrid, vector<double> &yvals,
   double startderiv, double endderiv) {
    Init (newGrid, yvals, startderiv, endderiv);
    Update();
  }
  /// Returns the value of the function at the ith grid point.
  inline double operator()(int i) const {
   return (y[i]);
  }
  /// Returns a reference to the value at the ith grid point.
  inline double & operator()(int i) {
   UpToDate = false;
   return (y[i]);
  }
  /// Trivial constructor
  CubicSpline() {
   UpToDate = false;
  }
};


inline double CubicSpline::operator()(double x) {
 if (!UpToDate)
  Update();

 GeneralGrid &X = grid;
#ifdef DEBUG
 if (x > X.End()) {
  if (x < (X.End() * 1.000000001))
   x = X.End();
  else {
   cerr << "x outside grid in CubicSpline.\n";
   cerr << "x = " << x << " X.End = " << X.End() << "\n";
   abort();
  }
 }
#endif
 int hi = X.ReverseMap(x)+1;
 int low = hi-1;
 if (low<0) {
  low = 0;
  hi = 1;
 }
 if (hi>(X.NumPoints()-1)) {
  hi = (X.NumPoints()-1);
  low = hi-1;
 }

 double h = X[hi] - X[low];
 double hinv = 1.0/h;
 double a = (X[hi]-x)*hinv;
 double b = (x-X[low])*hinv;
 double sixinv = 0.1666666666666666666;

 return (a*y[low] + b*y[hi] +
  ((a*a*a-a)*d2y[low]+(b*b*b-b)*d2y[hi])*(h*h*sixinv));
}


double CubicSpline::Deriv(double x) {
 if(!UpToDate)
  Update();

 GeneralGrid &X = grid;
 int hi = X.ReverseMap(x)+1;
 int low = hi-1;
 if (low<0) {
  low = 0;
  hi = 1;
 }
 if (hi>(X.NumPoints()-1)) {
  hi = (X.NumPoints()-1);
  low = hi-1;
 }

 double h = X[hi] - X[low];
 double hinv = 1.0/h;
 double a = (X[hi]-x)*hinv;
 double b = (x-X[low])*hinv;
 double sixinv = 0.1666666666666666666;

 return ((y[hi]-y[low])*hinv + (h*sixinv)*((3.0*b*b-1.0)*d2y[hi] -
  (3.0*a*a-1.0)*d2y[low]));
}


inline double CubicSpline::Deriv2(double x) {
 if(!UpToDate)
  Update();
 GeneralGrid &X = grid;
 int hi = X.ReverseMap(x)+1;
 int low = hi-1;
 if (low<0) {
  low = 0;
  hi = 1;
 }
 if (hi>(X.NumPoints()-1)) {
  hi = (X.NumPoints()-1);
  low = hi-1;
 }

 double h = X[hi] - X[low];
 double hinv = 1.0/h;
 double a = (X[hi]-x)*hinv;
 double b = (x-X[low])*hinv;

 return (a*d2y[low] + b*d2y[hi]);
}


inline double CubicSpline::Deriv3(double x) {
 if(!UpToDate)
   Update();
 GeneralGrid &X = grid;
 int hi = X.ReverseMap(x)+1;
 int low = hi-1;
 if (low<0) {
  low = 0;
  hi = 1;
 }
 if (hi>(X.NumPoints()-1)) {
  hi = (X.NumPoints()-1);
  low = hi-1;
 }

 double h = X[hi]-X[low];

 return ((d2y[hi]-d2y[low])/h);
}

#endif
