#ifndef INTERPOLATE_H
#define INTERPOLATE_H

double LinearInterpolate(double x, double xmin, double xmax, int nstops, const double* yarray)
{
  if (x<xmin) { x=xmin; }
  if (x>xmax) { x=xmax; }


  double ix = (nstops-1) * (x-xmin) / (xmax-xmin);
  //cout << "nstops=" << nstops << "\tx=" << x << "\txmin=" << xmin << "\txmax=" << xmax << "\tix = " << ix << endl;
  int i = int(ix);

  if (i==nstops-1) {
      return yarray[nstops-1];

    } // last element, no interp.
  double y0 = yarray[i];

  double y1 = yarray[i+1];
  return y0 + (y1-y0) * (ix-i);
}

double LinearInterpolate(double x, double xmin, double xmax, const vector<double> yvect) {
  return LinearInterpolate(x, xmin, xmax, yvect.size(), &yvect[0]);
}

#endif // INTERPOLATE_H

