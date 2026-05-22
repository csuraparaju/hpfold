#pragma once
#include <stdint.h>

/**
 * Lightweight rng state
 */
typedef struct {
    uint64_t state;
} rng_state;

/**
 * set the seed for the rng state. Defaults to an implementation
 * defined seed if seed is zero.
 */
void rng_seed(rng_state* rng, uint64_t seed);

/**
 * Returns a uniform random u64 integer
 */
uint64_t rng_u64(rng_state* rng);

/**
 * Returns a random uniform double in [0, 1]
 */
double rng_uniform(rng_state* rng);