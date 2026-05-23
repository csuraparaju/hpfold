#pragma once

#include "hp.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Interface for the host-side controller for the CUDA parallel tempering run
 */
typedef struct pt_sampler pt_sampler;

/**
 * Create a parallel tempering version of the metropolis hastings algorithm
 * with n_replicas, each running on a different temperature.
 * REQUIRES: \length(temperatures) == n_replicas
 * REQUIRES: is_sorted(temperatures)
 */
pt_sampler* pt_create(const char* seq, const double* temperatures, int n_replicas);

/**
 * Deallocate all memory associated with this parallel tempering data structure
 */
void pt_free(pt_sampler* s);

/**
 * Run n_swap rounds of the sampler. Each round has steps_per_swap MH steps and then
 * makes one swap attempt between each replica.
 */
void pt_run(pt_sampler* s, int steps_per_swap, int n_swaps);

/**
 * Return the best chain found so far
 */
const hp_chain* pt_best_chain(const pt_sampler* s);

/**
 * Return the acceptance rate of the swaps
 */
double pt_swap_acceptance_rate(const pt_sampler* s);

#ifdef __cplusplus
}
#endif