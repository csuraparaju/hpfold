#include <math.h>
#include <stdlib.h>
#include <stdbool.h>

#include "include/contracts.h"
#include "include/hp.h"
#include "include/mh.h"


// Forward declare the moves that this implementation uses.
bool hp_end_flip( hp_chain* chain, int aa_idx, rng_state* rng);
bool hp_corner_flip(hp_chain* chain, int aa_idx);
bool hp_crankshaft(hp_chain* chain, int aa_idx);

struct mh_sampler {
    hp_chain* chain;       // owned, current conformation
    hp_chain* best_chain;  // owned, lowest-energy conformation seen so far
    double temperature;
    rng_state rng;
    int steps_taken;
    int steps_accepted;
};

static bool is_mh_sampler(const mh_sampler* s) {
    if (s == NULL) return false;
    if (s->chain == NULL) return false;
    if (s->best_chain == NULL) return false;
    if (s->temperature <= 0.0) return false;
    if (s->steps_taken < 0) return false;
    if (s->steps_accepted < 0) return false;
    if (s->steps_accepted > s->steps_taken) return false;
    return true;
}

mh_sampler* mh_create(hp_chain* chain, double temperature) {
    REQUIRES(chain != NULL);
    REQUIRES(temperature > 0.0);

    mh_sampler* s = malloc(sizeof(mh_sampler));
    if (s == NULL) exit(EXIT_FAILURE);

    s->chain = chain;
    s->best_chain = hp_chain_clone(chain);
    s->temperature = temperature;
    s->steps_taken = 0;
    s->steps_accepted = 0;
    rng_seed(&s->rng, 0);
    ENSURES(is_mh_sampler(s));
    return s;
}

void mh_free(mh_sampler* s) {
    REQUIRES(is_mh_sampler(s));
    hp_chain_free(s->chain);
    hp_chain_free(s->best_chain);
    free(s);
}

void mh_step(mh_sampler* s) {
    REQUIRES(is_mh_sampler(s));

    hp_chain* chain = s->chain;
    int n = hp_chain_length(chain);

    // Pick move type uniformly
    // 0 = end_flip, 1 = corner_flip, 2 = crankshaft
    int move_type = (int)(rng_uniform(&s->rng) * 3.0);
    if (move_type > 2) move_type = 2;

    // Reject if chain is too short for this move type
    if (move_type == 1 && n < 3) { s->steps_taken++; return; }
    if (move_type == 2 && n < 4) { s->steps_taken++; return; }

    // Pick a valid residue for the chosen move type
    int aa_idx;
    if (move_type == 0) {
        aa_idx = (rng_uniform(&s->rng) < 0.5) ? 0 : n - 1;
    } else if (move_type == 1) {
        aa_idx = 1 + (int)(rng_uniform(&s->rng) * (n - 2));
        if (aa_idx >= n - 1) aa_idx = n - 2;
    } else {
        ASSERT(move_type == 2);
        aa_idx = 1 + (int)(rng_uniform(&s->rng) * (n - 3));
        if (aa_idx >= n - 2) aa_idx = n - 3;
    }

    // Save old position(s) and energy before the move so we can undo it
    // if MH rejects
    int old_energy = hp_chain_energy(chain);
    int old_x1, old_y1, old_x2, old_y2;
    hp_chain_get_coord(chain, aa_idx, &old_x1, &old_y1);
    if (move_type == 2) {
        hp_chain_get_coord(chain, aa_idx + 1, &old_x2, &old_y2);
    }

    // Propose and commit the move
    bool moved;
    if (move_type == 0) {
        moved = hp_end_flip(chain, aa_idx, &s->rng);
    }
    else if (move_type == 1) {
        moved = hp_corner_flip(chain, aa_idx);
    }
    else {
        moved = hp_crankshaft(chain, aa_idx);
    }

    s->steps_taken++;

    if (!moved) {
        ENSURES(is_mh_sampler(s));
        return;
    }

    // MH acceptance criterion
    int dE = hp_chain_energy(chain) - old_energy;
    bool accept = (dE <= 0) || (rng_uniform(&s->rng) < exp(-(double)dE / s->temperature));

    if (accept) {
        s->steps_accepted++;
        // Update best chain snapshot
        if (hp_chain_energy(chain) < hp_chain_energy(s->best_chain)) {
            hp_chain_free(s->best_chain);
            s->best_chain = hp_chain_clone(chain);
        }
    } else {
        // Reverse the move by committing it back to the old coordinates
        int new_x1, new_y1;
        hp_chain_get_coord(chain, aa_idx, &new_x1, &new_y1);

        if (move_type == 2) {
            int new_x2, new_y2;
            hp_chain_get_coord(chain, aa_idx + 1, &new_x2, &new_y2);
            hp_chain_commit_two_site_move(chain,
                aa_idx,     new_x1, new_y1, old_x1, old_y1,
                aa_idx + 1, new_x2, new_y2, old_x2, old_y2);
        } else {
            hp_chain_commit_move(chain, aa_idx, new_x1, new_y1, old_x1, old_y1);
        }
    }
    ENSURES(is_mh_sampler(s));
}

void mh_run(mh_sampler* s, int steps) {
    REQUIRES(is_mh_sampler(s));
    REQUIRES(steps >= 0);
    for (int i = 0; i < steps; i++) {
        mh_step(s);
    }
}

double mh_acceptance_rate(const mh_sampler* s) {
    REQUIRES(is_mh_sampler(s));
    if (s->steps_taken == 0) return 0.0;
    return (double)s->steps_accepted / s->steps_taken;
}

const hp_chain* mh_chain(const mh_sampler* s) {
    REQUIRES(is_mh_sampler(s));
    ENSURES(s->chain != NULL);
    return s->chain;
}

const hp_chain* mh_best_chain(const mh_sampler* s) {
    REQUIRES(is_mh_sampler(s));
    ENSURES(s->best_chain != NULL);
    return s->best_chain;
}
