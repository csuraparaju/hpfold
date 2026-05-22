#include <stdbool.h>
#include <stddef.h>
#include "rng.h"
#include "include/contracts.h"

/*
 * xorshift64 requires a non-zero state — zero is a fixed point that produces
 * zero forever.
 */
static bool is_rng_state(const rng_state* rng) {
    if (rng == NULL) return false;
    if (rng->state == 0) return false;
    return true;
}

void rng_seed(rng_state* rng, uint64_t seed) {
    REQUIRES(rng != NULL);
    rng->state = seed ? seed : 0xDEADBEEFCAFEBABEULL;
    ENSURES(is_rng_state(rng));
}

uint64_t rng_u64(rng_state* rng) {
    REQUIRES(is_rng_state(rng));

    uint64_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng->state = x;

    ENSURES(is_rng_state(rng));  // xorshift64 never maps non-zero to zero
    return x;
}

/*
 * Uniform double in [0, 1). Takes the top 53 bits of an xorshift64 output
 * and divides by 2^53, the largest power of two exactly representable in a double.
 */
double rng_uniform(rng_state* rng) {
    REQUIRES(is_rng_state(rng));
    return (rng_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}
