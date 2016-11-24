#ifndef HISTOGRAMS_HH
#define HISTOGRAMS_HH

#include <vector>

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
void save_histogram(const std::vector<double> &v, double dm, const char *filename);

#endif
