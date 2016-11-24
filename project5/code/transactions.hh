#ifndef TRANSACTIONS_HH
#define TRANSACTIONS_HH

#include <vector>

unsigned int seed_simulation();

/**
 * Simulate transactions.
 *
 * Arguments:
 *  size_t N    - Number of agents
 *  size_t K    - Number of transactions
 *  double m0   - Initial amount of money for each agent
 *  double l    - Saving fraction λ (default: 0.0)
 *  double a    - Neighbor preference parameter α (default: 0.0)
 *  double g    - Previous interaction preference parameter γ (default: 0.0)
 *
 * Returns:
 *  Vector with the amount of money each agent has after all transactions
**/
std::vector<double> simulate_transactions(size_t N, size_t K, double m0);
std::vector<double> simulate_transactions(size_t N, size_t K, double m0, double l);
std::vector<double> simulate_transactions(size_t N, size_t K, double m0, double l, double a);
std::vector<double> simulate_transactions(size_t N, size_t K, double m0, double l, double a, double g);

#endif
