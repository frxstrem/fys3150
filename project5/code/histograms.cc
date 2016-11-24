#include <iostream>
#include <cstdio>
#include <armadillo>
#include <random>

#include "histograms.hh"

using namespace std;
using namespace arma;

/**
 * Save histogram data to file.
 *
 * Will save file with four data points for each histogram bin:
 *  bin_start     Low edge of bin
 *  bin_end       High edge of bin
 *  count         Absolute count of values in bin
 *  relcount      Relative count (i.e., divided by total number of values)
 *
 * Arguments:
 *  v           - Vector containing all values
 *  dm          - Bin width
 *  filename    - Name of file to save to
**/
void save_histogram(const vector<double> &v, double dm, const char *filename) {
  // create armadillo vector from std::vector
  vec m(v);

  double mmin = 0, mmax = m.max();
  size_t nbins = ceil((mmax - mmin) / dm);
  size_t total = m.n_elem;

  vec bin_edges(nbins + 1);
  for(size_t i = 0; i <= nbins; i++)
    bin_edges[i] = mmin + dm * i;

  uvec bins = histc(m, bin_edges);

  // save histogram bins to file
  FILE *fp = fopen(filename, "w");
  fprintf(fp, "m_start\tm_end\tcount\trelcount\n");
  for(size_t i = 0; i < nbins; i++) {
    int count = bins(i);
    fprintf(fp, "%.3E\t%.3E\t%d\t%.3E\n", bin_edges(i), bin_edges(i + 1), count, (double)count / total);
  }
  fprintf(fp, "%.3E\t%.3E\t%d\t%.3E\n", bin_edges(nbins), bin_edges(nbins), 0, 0.0);
}
