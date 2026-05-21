#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "lib/contracts.h"
#include "hp.h"

#define MAX_N   64 // max chain length
#define LATTICE (2*MAX_N) // lattice is LATTICE x LATTICE, must be able to fit linear chain
#define ORIGIN  (LATTICE/2) // place amin acid 0 here so chain can grow in any direction

/**
 * Main data structure used to keep track of a Monte Carlo
 * trajectory for the HP folding model, as well as the helpers
 * and invariants on this data structure.
 */
struct mc_trajectory {
    int x[MAX_N]; // x-coordinate of each amino acid index
    int y[MAX_N]; // y-coordinate of each amino acid index
    int n; // actual chain length
    char aa_seq[MAX_N]; // 'H' for hydrophobic amino acid 'P' for polar amino acids
    int lattice_occ[LATTICE * LATTICE]; // occupancy lattice: -1 if point is empty, else amino acid index
    int energy; // current energy (negative count of H-H contacts)
};

#define LATTICE_IDX(x, y) ((x) * LATTICE + (y)) // row-major index into lattice_occ

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
        int j = traj->lattice_occ[LATTICE_IDX(nx, ny)];

        if (j == -1) continue;
        if (j == aa_idx) continue;
        if (j == aa_idx - 1 || j == aa_idx + 1) continue; // bonded neighbors are never contacts
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

    // Every non-empty site must hold a valid amino acid
    // index, and the total count must equal n
    int num_occupied = 0;
    for (int idx = 0; idx < LATTICE * LATTICE; idx++) {
        int b = traj->lattice_occ[idx];
        if (b == -1) continue;
        if (b < 0 || b >= traj->n) return false;
        num_occupied++;
    }
    if (num_occupied != traj->n) return false;

    // Check that each amino acid and its position is valid
    for (int i = 0; i < traj->n; i++) {
        if (traj->x[i] < 0 || traj->x[i] >= LATTICE) return false;
        if (traj->y[i] < 0 || traj->y[i] >= LATTICE) return false;
        if (traj->aa_seq[i] != 'H' && traj->aa_seq[i] != 'P') return false;

        int site = LATTICE_IDX(traj->x[i], traj->y[i]);
        if (traj->lattice_occ[site] != i) return false;
    }

    // Check that the chain is connected
    for (int i = 0; i < traj->n - 1; i++) {
        int mdist = abs(traj->x[i] - traj->x[i+1]) + abs(traj->y[i] - traj->y[i+1]);
        if (mdist != 1) return false;
    }

    if (traj->energy != compute_energy(traj)) return false;

    return true;
}

// Begin Interface Implementation

/**
 * Initialize traj with aa_seq laid out as a straight vertical chain
 * starting at ORIGIN, with all lattice sites cleared.
 */
mc_trajectory* mc_traj_create(const char* aa_seq) {
    REQUIRES(aa_seq != NULL);
    int n = (int)strlen(aa_seq);
    REQUIRES(n > 0 && n <= MAX_N);

    mc_trajectory* traj = malloc(sizeof(mc_trajectory));
    if (traj == NULL) exit(EXIT_FAILURE);

    traj->n = n;
    for (int i = 0; i < LATTICE * LATTICE; i++) traj->lattice_occ[i] = -1;

    for (int i = 0; i < n; i++) {
        traj->aa_seq[i] = aa_seq[i];
        traj->x[i] = ORIGIN;
        traj->y[i] = ORIGIN + i;
        traj->lattice_occ[LATTICE_IDX(ORIGIN, ORIGIN + i)] = i;
    }

    traj->energy = compute_energy(traj);
    ASSERT(traj->energy == 0); // sanity check, should be 0 for a linear chain

    ENSURES(is_mc_trajectory((traj)));
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

int mc_traj_run(mc_trajectory* traj, int steps, double temp) {
    REQUIRES(is_mc_trajectory(traj));
    REQUIRES(steps >= 0);
    REQUIRES(temp > 0.0);
    (void)temp;
    int best = traj->energy;
    for (int s = 0; s < steps; s++) {
        // TODO: attempt a random MC move and apply Metropolis acceptance
        if (traj->energy < best) best = traj->energy;
    }
    return best;
}

// End Interface Implementation

int main() {
    return 0;
}