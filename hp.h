#pragma once
#include <stdint.h>


typedef struct mc_trajectory mc_trajectory;

/**
 * Allocate and initialize a trajectory for the given null-terminated
 * HP sequence (e.g. "HPPHH").
 */
mc_trajectory* mc_traj_create(const char* seq);

/** Free all memory associated with traj. */
void mc_traj_free(mc_trajectory* traj);

/**
 * Run steps iterations of Metropolis Monte Carlo at the given temperature.
 * Returns the lowest energy observed over the run.
 */
int mc_traj_run(mc_trajectory* traj, int steps, double temp);

/** Current energy of the chain (negative H-H contact count). */
int mc_traj_energy(const mc_trajectory* traj);

/**
 * Write the x and y coordinates of each amino acid into xs[] and ys[].
 * Both arrays must be at least n elements, where n is the sequence length.
 */
void mc_traj_get_coords(const mc_trajectory* traj, int* xs, int* ys, int n);
