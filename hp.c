#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "lib/contracts.h"
#include "hp.h"

#define MAX_N   64  // max chain length
#define LATTICE (2*MAX_N) // lattice is LATTICE x LATTICE, must fit a linear chain
#define ORIGIN  (LATTICE/2) // place amino acid 0 here so chain can grow in any direction

#define LATTICE_IDX(x, y) ((x) * LATTICE + (y)) // row-major index into lattice_occ

// Default RNG seed used when none is supplied by the interface.
// xorshift64 cannot accept a zero state.
#define DEFAULT_SEED 0xDEADBEEFCAFEBABEULL

/**
 * Main data structure used to keep track of a Monte Carlo trajectory for the
 * HP folding model, as well as the helpers and invariants on this data structure.
 */
struct mc_trajectory {
    int x[MAX_N];                       // x-coordinate of each amino acid index
    int y[MAX_N];                       // y-coordinate of each amino acid index
    int n;                              // actual chain length
    char aa_seq[MAX_N];                 // 'H' for hydrophobic, 'P' for polar
    int lattice_occ[LATTICE * LATTICE]; // amino acid index if it is at the specific (x, y) coord, else -1
    int energy;                         // current energy (negative count of H-H contacts)
    uint64_t rng_state;                 // state for the xorshift64 RNG
};

/**
 * Count non-bonded H-H contacts amino acid aa_idx would have if placed at (cx, cy).
 * Used to compute energy changes incrementally without a full rescan.
 */
static int count_contacts_at(const mc_trajectory* traj, int aa_idx, int cx, int cy) {
    REQUIRES(traj != NULL);
    REQUIRES(aa_idx >= 0 && aa_idx < traj->n);
    REQUIRES(cx >= 0 && cx < LATTICE && cy >= 0 && cy < LATTICE);

    static const int dx[4] = {1, -1, 0, 0};
    static const int dy[4] = {0, 0, 1, -1};
    int contacts = 0;

    if (traj->aa_seq[aa_idx] != 'H') return 0;

    for (int d = 0; d < 4; d++) {
        int nx = cx + dx[d];
        int ny = cy + dy[d];
        if (nx < 0 || nx >= LATTICE || ny < 0 || ny >= LATTICE) continue;  // off-lattice
        int j = traj->lattice_occ[LATTICE_IDX(nx, ny)];

        if (j == -1) continue;
        if (j == aa_idx) continue;
        if (j == aa_idx - 1 || j == aa_idx + 1) continue;  // bonded neighbors are never contacts
        if (traj->aa_seq[j] != 'H') continue;

        contacts++;
    }
    return contacts;
}

/**
 * Sum all H-H contacts in the chain. Divides by 2 because each contact
 * is counted once from each endpoint.
 */
static int compute_energy(const mc_trajectory* traj) {
    REQUIRES(traj != NULL);
    REQUIRES(traj->n > 0 && traj->n <= MAX_N);
    int contacts = 0;

    for (int i = 0; i < traj->n; i++) {
        contacts += count_contacts_at(traj, i, traj->x[i], traj->y[i]);
    }
    return -contacts / 2;
}

/**
 * Main data structure invariant on mc_trajectory.
 */
static bool is_mc_trajectory(const mc_trajectory* traj) {
    if (traj == NULL) return false;
    if (traj->n <= 0 || traj->n > MAX_N) return false;
    if (traj->rng_state == 0) return false;  // xorshift64 invariant

    // Every non-empty site must hold a valid amino acid index,
    // and the total count must equal n.
    int num_occupied = 0;
    for (int idx = 0; idx < LATTICE * LATTICE; idx++) {
        int b = traj->lattice_occ[idx];
        if (b == -1) continue;
        if (b < 0 || b >= traj->n) return false;
        num_occupied++;
    }
    if (num_occupied != traj->n) return false;

    // Check that each amino acid and its position is valid.
    for (int i = 0; i < traj->n; i++) {
        if (traj->x[i] < 0 || traj->x[i] >= LATTICE) return false;
        if (traj->y[i] < 0 || traj->y[i] >= LATTICE) return false;
        if (traj->aa_seq[i] != 'H' && traj->aa_seq[i] != 'P') return false;

        int site = LATTICE_IDX(traj->x[i], traj->y[i]);
        if (traj->lattice_occ[site] != i) return false;
    }

    // Check that the chain is connected.
    for (int i = 0; i < traj->n - 1; i++) {
        int mdist = abs(traj->x[i] - traj->x[i+1]) + abs(traj->y[i] - traj->y[i+1]);
        if (mdist != 1) return false;
    }

    if (traj->energy != compute_energy(traj)) return false;

    return true;
}

// ============================================================
// Interface implementation
// ============================================================

/**
 * xorshift64 PRNG. State must never be zero
 */
static uint64_t xorshift64(uint64_t* state) {
    REQUIRES(state != NULL);
    REQUIRES(*state != 0);

    uint64_t x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    return *state = x;
}

/**
 * Uniform double in [0, 1). Takes the top 53 bits of an xorshift64 output
 * and divides by 2^53, the largest power of two exactly representable in a double.
 */
static double rand_uniform(uint64_t* state) {
    REQUIRES(state != NULL);
    return (xorshift64(state) >> 11) * (1.0 / 9007199254740992.0);
}

/**
 * Allocate and initialize a trajectory for aa_seq, laid out as a straight
 * line along the y-axis starting at (ORIGIN, ORIGIN).
 */
mc_trajectory* mc_traj_create(const char* aa_seq) {
    REQUIRES(aa_seq != NULL);
    int n = (int)strlen(aa_seq);
    REQUIRES(n > 0 && n <= MAX_N);

    mc_trajectory* traj = malloc(sizeof(mc_trajectory));
    if (traj == NULL) exit(EXIT_FAILURE);

    traj->n = n;
    traj->rng_state = DEFAULT_SEED;

    for (int i = 0; i < LATTICE * LATTICE; i++) traj->lattice_occ[i] = -1;

    for (int i = 0; i < n; i++) {
        traj->aa_seq[i] = aa_seq[i];
        traj->x[i] = ORIGIN;
        traj->y[i] = ORIGIN + i;
        traj->lattice_occ[LATTICE_IDX(ORIGIN, ORIGIN + i)] = i;
    }

    traj->energy = compute_energy(traj);
    ASSERT(traj->energy == 0);  // straight chain has no H-H contacts

    ENSURES(is_mc_trajectory(traj));
    return traj;
}

void mc_traj_free(mc_trajectory* traj) {
    REQUIRES(is_mc_trajectory(traj));
    free(traj);
}

int mc_traj_energy(const mc_trajectory* traj) {
    REQUIRES(is_mc_trajectory(traj));
    return traj->energy;
}

void mc_traj_get_coords(const mc_trajectory* traj, int* xs, int* ys, int n) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(xs != NULL && ys != NULL);
    REQUIRES(n == traj->n);
    for (int i = 0; i < n; i++) {
        xs[i] = traj->x[i];
        ys[i] = traj->y[i];
    }
}


/**
 * Attempt an end-flip on the specified endpoint (amino acid 0 or amino acid n-1).
 * Picks one of the 3 candidate positions uniformly, runs the
 * Metropolis accept/reject pipeline, and commits if accepted.
 * Returns true if the move was accepted, false if rejected (for any reason:
 * off-lattice, self-avoidance violation, or Metropolis rejection).
 */
static bool make_end_flip_move(mc_trajectory* traj, int aa_idx, double T) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(aa_idx == 0 || aa_idx == traj->n - 1); // End amino acid
    REQUIRES(T > 0.0);

    // Find the sequence-neighbor the endpoint pivots around.
    int nbor_idx = (aa_idx == 0) ? 1 : aa_idx - 1;

    int aa_x   = traj->x[aa_idx];
    int aa_y   = traj->y[aa_idx];
    int nbor_x = traj->x[nbor_idx];
    int nbor_y = traj->y[nbor_idx];

    ASSERT(abs(aa_x - nbor_x) + abs(aa_y - nbor_y) == 1);

    static const int dx[4] = {1, -1, 0, 0};
    static const int dy[4] = {0, 0, 1, -1};
    int cand_x[3], cand_y[3];
    int n_cand = 0;
    for (int d = 0; d < 4; d++) {
        int cx = nbor_x + dx[d];
        int cy = nbor_y + dy[d];
        if (cx == aa_x && cy == aa_y) continue;
        cand_x[n_cand] = cx;
        cand_y[n_cand] = cy;
        n_cand++;
    }
    ASSERT(n_cand == 3);

    int choice = (int)(rand_uniform(&traj->rng_state) * 3);
    int new_x = cand_x[choice];
    int new_y = cand_y[choice];

    // Move is out of bounds, ignore
    if (new_x < 0 || new_x >= LATTICE || new_y < 0 || new_y >= LATTICE) {
        return false;
    }

    // Move will over ride existing amino acid, ignore
    if (traj->lattice_occ[LATTICE_IDX(new_x, new_y)] != -1) {
        return false;
    }

    // Energy delta. count_contacts_at excludes bonded neighbors and self,
    // so we can safely compute new_contacts even though aa_idx is still
    // recorded at (aa_x, aa_y) in lattice_occ.
    int old_contacts = count_contacts_at(traj, aa_idx, aa_x, aa_y);
    int new_contacts = count_contacts_at(traj, aa_idx, new_x, new_y);
    int dE = -(new_contacts - old_contacts);

    // Metropolis acceptance condition
    if (dE > 0 && rand_uniform(&traj->rng_state) >= exp(-dE / T)) {
        return false;
    }

    // Accept this move and update accordingly
    traj->lattice_occ[LATTICE_IDX(aa_x, aa_y)] = -1;
    traj->lattice_occ[LATTICE_IDX(new_x, new_y)] = aa_idx;
    traj->x[aa_idx] = new_x;
    traj->y[aa_idx] = new_y;
    traj->energy += dE;

    ENSURES(is_mc_trajectory(traj));
    return true;
}

/**
 * Attempt a corner-flip on the specified interior amino acid.
 * Eligible iff the amino acid sits at one corner of a unit square formed
 * with its two sequence-neighbors; the move flips it to the opposite corner.
 * Returns true if accepted, false if rejected (ineligible, off-lattice,
 * self-avoidance violation, or Metropolis rejection).
 */
static bool make_corner_flip_move(mc_trajectory* traj, int aa_idx, double T) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(aa_idx > 0 && aa_idx < traj->n - 1);  // interior amino acid
    REQUIRES(T > 0.0);

    int prev_x = traj->x[aa_idx - 1];
    int prev_y = traj->y[aa_idx - 1];
    int next_x = traj->x[aa_idx + 1];
    int next_y = traj->y[aa_idx + 1];

    // corner flip move doesn't apple because aa_idx is in a straight
    // segment, ignore
    if (prev_x == next_x || prev_y == next_y) return false;

    int aa_x = traj->x[aa_idx];
    int aa_y = traj->y[aa_idx];

    // aa_idx is currently at one of (prev_x, next_y) or
    // (next_x, prev_y), the candidate is the other
    int new_x = (aa_x == prev_x) ? next_x : prev_x;
    int new_y = (aa_y == prev_y) ? next_y : prev_y;

    ASSERT(abs(new_x - prev_x) + abs(new_y - prev_y) == 1);
    ASSERT(abs(new_x - next_x) + abs(new_y - next_y) == 1);
    ASSERT(new_x != aa_x || new_y != aa_y);
    ASSERT(new_x >= 0 && new_x < LATTICE && new_y >= 0 && new_y < LATTICE);

    // Move will over ride existing amino acid, ignore
    if (traj->lattice_occ[LATTICE_IDX(new_x, new_y)] != -1) return false;

    // Energy delta computation
    int old_contacts = count_contacts_at(traj, aa_idx, aa_x, aa_y);
    int new_contacts = count_contacts_at(traj, aa_idx, new_x, new_y);
    int dE = -(new_contacts - old_contacts);

    // Metropolis acceptance condition
    if (dE > 0 && rand_uniform(&traj->rng_state) >= exp(-dE / T)) return false;

    // Accept this move and update accordingly
    traj->lattice_occ[LATTICE_IDX(aa_x, aa_y)] = -1;
    traj->lattice_occ[LATTICE_IDX(new_x, new_y)] = aa_idx;
    traj->x[aa_idx] = new_x;
    traj->y[aa_idx] = new_y;
    traj->energy += dE;

    ENSURES(is_mc_trajectory(traj));
    return true;
}

/**
 * Attempt a crankshaft on the pair (aa_idx, aa_idx+1).
 * Eligible iff amino acids aa_idx-1 and aa_idx+2 are lattice-neighbors,
 * forming a U-shape that amino acids aa_idx and aa_idx+1 hang off of.
 * Reflects both amino acids across the line connecting aa_idx-1 and aa_idx+2.
 * Returns true if accepted, false if rejected.
 */
static bool make_crankshaft_move(mc_trajectory* traj, int aa_idx, double T) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(aa_idx >= 1 && aa_idx <= traj->n - 3);  // need i-1 and i+2 in bounds
    REQUIRES(T > 0.0);

    int ax = traj->x[aa_idx - 1], ay = traj->y[aa_idx - 1];   // anchor 1
    int bx = traj->x[aa_idx + 2], by = traj->y[aa_idx + 2];   // anchor 2

    // anchors are not lattice points, ignore
    int adist = abs(ax - bx) + abs(ay - by);
    if (adist != 1) return false;

    int i_x  = traj->x[aa_idx];
    int i_y  = traj->y[aa_idx];
    int i1_x = traj->x[aa_idx + 1];
    int i1_y = traj->y[aa_idx + 1];

    // Reflect the anchors
    int new_i_x, new_i_y, new_i1_x, new_i1_y;
    if (ax == bx) {
        // Anchors share x, reflect across the vertical line x = ax.
        new_i_x  = 2 * ax - i_x;
        new_i_y  = i_y;

        new_i1_x = 2 * ax - i1_x;
        new_i1_y = i1_y;
    } else {
        // Anchors share y, reflect across horizontal line y = ay.
        ASSERT(ay == by);
        new_i_x  = i_x;
        new_i_y  = 2 * ay - i_y;

        new_i1_x = i1_x;
        new_i1_y = 2 * ay - i1_y;
    }

    ASSERT(new_i_x  >= 0 && new_i_x  < LATTICE && new_i_y  >= 0 && new_i_y  < LATTICE);
    ASSERT(new_i1_x >= 0 && new_i1_x < LATTICE && new_i1_y >= 0 && new_i1_y < LATTICE);

    // if the moving amino acids are already on the line of reflection,
    // they don't move, ignore
    if (new_i_x == i_x && new_i_y == i_y) return false;

    // Move will over ride existing amino acid, ignore
    int occ_new_i  = traj->lattice_occ[LATTICE_IDX(new_i_x,  new_i_y)];
    int occ_new_i1 = traj->lattice_occ[LATTICE_IDX(new_i1_x, new_i1_y)];
    if (occ_new_i  != -1 && occ_new_i  != aa_idx && occ_new_i  != aa_idx + 1) return false;
    if (occ_new_i1 != -1 && occ_new_i1 != aa_idx && occ_new_i1 != aa_idx + 1) return false;

    // Energy delta. Note: count_contacts_at(i) excludes the bonded neighbor
    // i+1 (and vice versa), so the i--(i+1) contacts are not counted
    int old_contacts = count_contacts_at(traj, aa_idx, i_x,  i_y)
                     + count_contacts_at(traj, aa_idx + 1, i1_x, i1_y);

    // Tentatively apply the move so we can compute the contacts
    traj->lattice_occ[LATTICE_IDX(i_x,  i_y)]  = -1;
    traj->lattice_occ[LATTICE_IDX(i1_x, i1_y)] = -1;
    traj->lattice_occ[LATTICE_IDX(new_i_x,  new_i_y)]  = aa_idx;
    traj->lattice_occ[LATTICE_IDX(new_i1_x, new_i1_y)] = aa_idx + 1;
    traj->x[aa_idx] = new_i_x;
    traj->y[aa_idx] = new_i_y;
    traj->x[aa_idx + 1] = new_i1_x;
    traj->y[aa_idx + 1] = new_i1_y;

    int new_contacts = count_contacts_at(traj, aa_idx,     new_i_x,  new_i_y)
                     + count_contacts_at(traj, aa_idx + 1, new_i1_x, new_i1_y);
    int dE = -(new_contacts - old_contacts);

    // Metropolis acceptance condition
    if (dE > 0 && rand_uniform(&traj->rng_state) >= exp(-dE / T)) {
        // Reject because the probability
        traj->lattice_occ[LATTICE_IDX(new_i_x,  new_i_y)]  = -1;
        traj->lattice_occ[LATTICE_IDX(new_i1_x, new_i1_y)] = -1;
        traj->lattice_occ[LATTICE_IDX(i_x,  i_y)]  = aa_idx;
        traj->lattice_occ[LATTICE_IDX(i1_x, i1_y)] = aa_idx + 1;
        traj->x[aa_idx]     = i_x;   traj->y[aa_idx]     = i_y;
        traj->x[aa_idx + 1] = i1_x;  traj->y[aa_idx + 1] = i1_y;
        return false;
    }

    // Accept and return true
    traj->energy += dE;
    ENSURES(is_mc_trajectory(traj));
    return true;
}

int mc_traj_run(mc_trajectory* traj, int steps, double temp) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(traj->n >= 4); // Moves aren't defined for smaller aa chains
    REQUIRES(steps >= 0);
    REQUIRES(temp > 0.0);

    int best = traj->energy;

    int best_x[MAX_N];
    int best_y[MAX_N];
    for (int i = 0; i < traj->n; i++) {
        best_x[i] = traj->x[i];
        best_y[i] = traj->y[i];
    }

    int n_accepted = 0;

    for (int s = 0; s < steps; s++) {
        double r = rand_uniform(&traj->rng_state);
        if (r < 1.0/3.0 ) {
            int end = (rand_uniform(&traj->rng_state) < 0.5) ? 0 : traj->n - 1;
            if (make_end_flip_move(traj, end, temp)) n_accepted++;
        } else if (r < 2.0/3.0) {
            int idx = 1 + (int)(rand_uniform(&traj->rng_state) * (traj->n - 2));
            if (make_corner_flip_move(traj, idx, temp)) n_accepted++;
        } else {
            int idx = 1 + (int)(rand_uniform(&traj->rng_state) * (traj->n - 3));
            if (make_crankshaft_move(traj, idx, temp)) n_accepted++;
        }
        if (traj->energy < best) {
            best = traj->energy;
            for (int i = 0; i < traj->n; i++) {
                best_x[i] = traj->x[i];
                best_y[i] = traj->y[i];
            }
        }

        if (s % 1000 == 0) {
            ASSERT(traj->energy == compute_energy(traj));
        }

    }
    fprintf(stdout, "accepted %d/%d (%.1f%%) min=%d\n", n_accepted, steps, 100.0*n_accepted/steps, best);

    // Restore to the best configuration
    for (int i = 0; i < LATTICE*LATTICE; i++)
        traj->lattice_occ[i] = -1;

    for (int i = 0; i < traj->n; i++) {
        traj->x[i] = best_x[i];
        traj->y[i] = best_y[i];

        traj->lattice_occ[LATTICE_IDX(best_x[i], best_y[i])] = i;
    }

    traj->energy = best;

    ENSURES(is_mc_trajectory(traj));
    return best;
}

void mc_traj_print(const mc_trajectory* traj) {
    REQUIRES(is_mc_trajectory(traj));

    // Expanded canvas:
    // amino acids at (2x, 4y) so labels fit nicely
    const int ROWS = 2 * LATTICE - 1;
    const int COLS = 4 * LATTICE - 1;

    char canvas[ROWS][COLS + 1];

    // Fill with spaces
    for (int r = 0; r < ROWS; r++) {
        for (int c = 0; c < COLS; c++) {
            canvas[r][c] = ' ';
        }
        canvas[r][COLS] = '\0';
    }

    // Draw amino acids and bonds
    for (int i = 0; i < traj->n; i++) {
        int x = traj->x[i];
        int y = traj->y[i];

        // Put lattice point on expanded grid
        int row = 2 * x;
        int col = 4 * y;

        // Print amino acid label (e.g. H3)
        char label[8];
        snprintf(label, sizeof(label), "%c", traj->aa_seq[i]);

        for (int k = 0; label[k] != '\0' && col+k < COLS; k++) {
            canvas[row][col+k] = label[k];
        }

        // Draw bond to next amino acid
        if (i < traj->n - 1) {
            int nx = traj->x[i+1];
            int ny = traj->y[i+1];

            int drow = 2*(nx - x);
            int dcol = 4*(ny - y);

            // Vertical bond
            if (drow != 0) {
                canvas[row + drow/2][col] = '|';
            }

            // Horizontal bond
            if (dcol != 0) {
                int step = (dcol > 0) ? 1 : -1;
                for (int c = col + step;
                     c != col + dcol;
                     c += step) {
                    canvas[row][c] = '-';
                }
            }
        }
    }

    // Find bounding box so we don't print huge whitespace
    int min_r = ROWS, max_r = 0;
    int min_c = COLS, max_c = 0;

    for (int r = 0; r < ROWS; r++) {
        for (int c = 0; c < COLS; c++) {
            if (canvas[r][c] != ' ') {
                if (r < min_r) min_r = r;
                if (r > max_r) max_r = r;
                if (c < min_c) min_c = c;
                if (c > max_c) max_c = c;
            }
        }
    }

    printf("\nEnergy = %d\n", traj->energy);

    for (int r = min_r; r <= max_r; r++) {
        for (int c = min_c; c <= max_c; c++) {
            putchar(canvas[r][c]);
        }
        putchar('\n');
    }
    putchar('\n');
}
