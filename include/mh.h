#pragma once

#include "hp.h"
#include "rng.h"

typedef struct mh_sampler mh_sampler;

/**
 * Create a metropolis hastings algorithm sampler for the given
 * chain and the temperature
 */
mh_sampler* mh_create(hp_chain* chain, double temperature);

/**
 * Free up all resources for this sampler
 */
void mh_free(mh_sampler* sampler);

/**
 * Perform one step of the metropolis hastings algorithm
 */
void mh_step(mh_sampler* sampler);

/**
 * Run steps number of steps for the algorithm
 */
void mh_run(mh_sampler* sampler, int steps);

/**
 * Get the acceptance rate of all the moves
 */
double mh_acceptance_rate(const mh_sampler* sampler);

/**
 * Return a read-only pointer to the sampler's current chain.
 */
const hp_chain* mh_chain(const mh_sampler* sampler);

/**
 * Return a read-only pointer to the lowest-energy chain seen so far.
 */
const hp_chain* mh_best_chain(const mh_sampler* sampler);