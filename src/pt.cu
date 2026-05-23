#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/hp.h"
#include "include/pt.h"
#include "include/rng.h"

/*
 * Chain state on device. Duplicated from src/hp.c
 */
struct hp_chain {
    int x[HP_MAX_N];
    int y[HP_MAX_N];
    int n;
    char aa_seq[HP_MAX_N];
    int energy;
};

__device__ int d_residue_at(const hp_chain* chain, int x, int y) {
    for (int i = 0; i < chain->n; i++) {
        if (chain->x[i] == x && chain->y[i] == y) return i;
    }
    return -1;
}

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,      \
                    cudaGetErrorString(err));                                  \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)

/*
 * Advance the per-replica xorshift64 state and return the next u64 value.
 */
__device__ uint64_t d_rng_u64(rng_state* rng) {
    uint64_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng->state = x;
    return x;
}

/*
 * Draw a uniform double in [0, 1) on the device.
 */
__device__ double d_rng_uniform(rng_state* rng) {
    return (double)(d_rng_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}

/*
 * Count non-bonded H-H contacts amino acid aa_idx would have at lattice site (cx, cy).
 */
__device__ int d_hp_chain_contacts_at(const hp_chain* chain, int aa_idx, int cx, int cy) {
    const int dx[4] = {1, -1, 0, 0};
    const int dy[4] = {0, 0, 1, -1};
    int contacts = 0;

    if (chain->aa_seq[aa_idx] != 'H') return 0;

    for (int d = 0; d < 4; d++) {
        int nx = cx + dx[d];
        int ny = cy + dy[d];
        if (nx < 0 || nx >= HP_LATTICE || ny < 0 || ny >= HP_LATTICE) continue;
        int j = d_residue_at(chain, nx, ny);
        if (j == -1) continue;
        if (j == aa_idx) continue;
        if (j == aa_idx - 1 || j == aa_idx + 1) continue;
        if (chain->aa_seq[j] != 'H') continue;
        contacts++;
    }
    return contacts;
}

/*
 * Read the (x, y) lattice coordinates of amino acid i.
 */
__device__ void d_hp_chain_get_coord(const hp_chain* chain, int i, int* x, int* y) {
    *x = chain->x[i];
    *y = chain->y[i];
}

/*
 * Return which amino acid occupies (x, y), or -1 if the site is empty.
 */
__device__ int d_hp_chain_site_occ(const hp_chain* chain, int x, int y) {
    return d_residue_at(chain, x, y);
}

/*
 * Move one amino acid on the lattice and update coordinates and energy.
 */
__device__ void d_hp_chain_commit_move(hp_chain* chain, int aa_idx,
                                       int old_x, int old_y, int new_x, int new_y) {
    int old_contacts = d_hp_chain_contacts_at(chain, aa_idx, old_x, old_y);
    int new_contacts = d_hp_chain_contacts_at(chain, aa_idx, new_x, new_y);
    int dE = -(new_contacts - old_contacts);

    chain->x[aa_idx] = new_x;
    chain->y[aa_idx] = new_y;
    chain->energy += dE;
}

/*
 * Move two amino acids atomically and apply the combined energy change.
 */
__device__ void d_hp_chain_commit_two_site_move(hp_chain* chain,
    int aa_idx1, int old_x1, int old_y1, int new_x1, int new_y1,
    int aa_idx2, int old_x2, int old_y2, int new_x2, int new_y2) {
    int old_contacts = d_hp_chain_contacts_at(chain, aa_idx1, old_x1, old_y1) +
                       d_hp_chain_contacts_at(chain, aa_idx2, old_x2, old_y2);

    chain->x[aa_idx1] = new_x1;
    chain->y[aa_idx1] = new_y1;
    chain->x[aa_idx2] = new_x2;
    chain->y[aa_idx2] = new_y2;

    int new_contacts = d_hp_chain_contacts_at(chain, aa_idx1, new_x1, new_y1) +
                       d_hp_chain_contacts_at(chain, aa_idx2, new_x2, new_y2);

    chain->energy += -(new_contacts - old_contacts);
}

/*
 * True if (x, y) lies inside the HP lattice bounds.
 */
__device__ bool d_in_bounds(int x, int y) {
    return x >= 0 && x < HP_LATTICE && y >= 0 && y < HP_LATTICE;
}

/*
 * Propose and apply an end-flip move on endpoint aa_idx; uses rng to pick a candidate site.
 */
__device__ bool d_hp_end_flip(hp_chain* chain, int aa_idx, rng_state* rng) {
    const int DX[4] = {1, -1, 0, 0};
    const int DY[4] = {0, 0, 1, -1};

    int nbor_idx = (aa_idx == 0) ? 1 : aa_idx - 1;
    int aa_x, aa_y, nbor_x, nbor_y;
    d_hp_chain_get_coord(chain, aa_idx, &aa_x, &aa_y);
    d_hp_chain_get_coord(chain, nbor_idx, &nbor_x, &nbor_y);

    int cand_x[3];
    int cand_y[3];
    int n_cand = 0;

    for (int d = 0; d < 4; d++) {
        int cx = nbor_x + DX[d];
        int cy = nbor_y + DY[d];
        if (cx == aa_x && cy == aa_y) continue;
        cand_x[n_cand] = cx;
        cand_y[n_cand] = cy;
        n_cand++;
    }

    int choice = (int)(d_rng_uniform(rng) * 3.0);
    if (choice < 0) choice = 0;
    if (choice > 2) choice = 2;

    int new_x = cand_x[choice];
    int new_y = cand_y[choice];

    if (!d_in_bounds(new_x, new_y)) return false;
    if (d_hp_chain_site_occ(chain, new_x, new_y) != -1) return false;

    d_hp_chain_commit_move(chain, aa_idx, aa_x, aa_y, new_x, new_y);
    return true;
}

/*
 * Propose and apply a corner-flip move on interior amino acid aa_idx.
 */
__device__ bool d_hp_corner_flip(hp_chain* chain, int aa_idx) {
    int prev_x, prev_y, next_x, next_y;
    d_hp_chain_get_coord(chain, aa_idx - 1, &prev_x, &prev_y);
    d_hp_chain_get_coord(chain, aa_idx + 1, &next_x, &next_y);

    if (prev_x == next_x || prev_y == next_y) return false;

    int aa_x, aa_y;
    d_hp_chain_get_coord(chain, aa_idx, &aa_x, &aa_y);

    int new_x = (aa_x == prev_x) ? next_x : prev_x;
    int new_y = (aa_y == prev_y) ? next_y : prev_y;

    if (new_x == aa_x && new_y == aa_y) return false;
    if (d_hp_chain_site_occ(chain, new_x, new_y) != -1) return false;

    d_hp_chain_commit_move(chain, aa_idx, aa_x, aa_y, new_x, new_y);
    return true;
}

/*
 * Propose and apply a crankshaft move on the consecutive pair (aa_idx, aa_idx + 1).
 */
__device__ bool d_hp_crankshaft(hp_chain* chain, int aa_idx) {
    int ax, ay, bx, by;
    d_hp_chain_get_coord(chain, aa_idx - 1, &ax, &ay);
    d_hp_chain_get_coord(chain, aa_idx + 2, &bx, &by);

    if (abs(ax - bx) + abs(ay - by) != 1) return false;

    int i_x, i_y, i1_x, i1_y;
    d_hp_chain_get_coord(chain, aa_idx, &i_x, &i_y);
    d_hp_chain_get_coord(chain, aa_idx + 1, &i1_x, &i1_y);

    int new_i_x, new_i_y, new_i1_x, new_i1_y;

    if (ax == bx) {
        new_i_x = 2 * ax - i_x;
        new_i_y = i_y;
        new_i1_x = 2 * ax - i1_x;
        new_i1_y = i1_y;
    } else {
        new_i_x = i_x;
        new_i_y = 2 * ay - i_y;
        new_i1_x = i1_x;
        new_i1_y = 2 * ay - i1_y;
    }

    if (new_i_x == i_x && new_i_y == i_y && new_i1_x == i1_x && new_i1_y == i1_y)
        return false;

    if (!d_in_bounds(new_i_x, new_i_y) || !d_in_bounds(new_i1_x, new_i1_y)) return false;
    if (d_hp_chain_site_occ(chain, new_i_x, new_i_y) != -1) return false;
    if (d_hp_chain_site_occ(chain, new_i1_x, new_i1_y) != -1) return false;

    d_hp_chain_commit_two_site_move(chain,
        aa_idx, i_x, i_y, new_i_x, new_i_y,
        aa_idx + 1, i1_x, i1_y, new_i1_x, new_i1_y);
    return true;
}

/*
 * One Metropolis–Hastings step: pick a move, apply it, accept or undo by temperature.
 */
__device__ void d_mh_step(hp_chain* chain, rng_state* rng, double temperature) {
    int n = chain->n;

    int move_type = (int)(d_rng_uniform(rng) * 3.0);
    if (move_type > 2) move_type = 2;

    if (move_type == 1 && n < 3) return;
    if (move_type == 2 && n < 4) return;

    int aa_idx;
    if (move_type == 0) {
        aa_idx = (d_rng_uniform(rng) < 0.5) ? 0 : n - 1;
    } else if (move_type == 1) {
        aa_idx = 1 + (int)(d_rng_uniform(rng) * (n - 2));
        if (aa_idx >= n - 1) aa_idx = n - 2;
    } else {
        aa_idx = 1 + (int)(d_rng_uniform(rng) * (n - 3));
        if (aa_idx >= n - 2) aa_idx = n - 3;
    }

    int old_energy = chain->energy;
    int old_x1, old_y1, old_x2, old_y2;
    d_hp_chain_get_coord(chain, aa_idx, &old_x1, &old_y1);
    if (move_type == 2) {
        d_hp_chain_get_coord(chain, aa_idx + 1, &old_x2, &old_y2);
    }

    bool moved;
    if (move_type == 0) {
        moved = d_hp_end_flip(chain, aa_idx, rng);
    } else if (move_type == 1) {
        moved = d_hp_corner_flip(chain, aa_idx);
    } else {
        moved = d_hp_crankshaft(chain, aa_idx);
    }

    if (!moved) return;

    int dE = chain->energy - old_energy;
    bool accept = (dE <= 0) || (d_rng_uniform(rng) < exp(-(double)dE / temperature));

    if (!accept) {
        int new_x1, new_y1;
        d_hp_chain_get_coord(chain, aa_idx, &new_x1, &new_y1);
        if (move_type == 2) {
            int new_x2, new_y2;
            d_hp_chain_get_coord(chain, aa_idx + 1, &new_x2, &new_y2);
            d_hp_chain_commit_two_site_move(chain,
                aa_idx, new_x1, new_y1, old_x1, old_y1,
                aa_idx + 1, new_x2, new_y2, old_x2, old_y2);
        } else {
            d_hp_chain_commit_move(chain, aa_idx, new_x1, new_y1, old_x1, old_y1);
        }
    }
}

/*
 * Run steps_per_swap MH steps on each replica; one CUDA thread per replica.
 */
__global__ void pt_mh_kernel(hp_chain* d_chains, double* d_temps,
                               rng_state* d_rngs, int n, int steps) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;

    hp_chain chain = d_chains[r];
    rng_state rng = d_rngs[r];
    double temp = d_temps[r];

    for (int s = 0; s < steps; s++) {
        d_mh_step(&chain, &rng, temp);
    }

    d_chains[r] = chain;
    d_rngs[r] = rng;
}

/*
 * Write each replica's energy into d_energies for host-side temperature swaps.
 */
__global__ void pt_energies_kernel(const hp_chain* d_chains, int* d_energies, int n) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n) return;
    d_energies[r] = d_chains[r].energy;
}

/*
 * Host-side state for a parallel tempering run: device replicas, host swap buffers, best chain.
 */
struct pt_sampler {
    hp_chain* d_chains;
    double* d_temps;
    rng_state* d_rngs;
    int* d_energies;

    int* h_energies;
    double* h_temps;

    hp_chain h_best;
    int n_replicas;
    long swap_attempts;
    long swap_accepted;

    rng_state swap_rng;
};

/*
 * Number of thread blocks needed to launch one thread per replica.
 */
static int pt_launch_blocks(int n) {
    const int threads = 256;
    return (n + threads - 1) / threads;
}

/*
 * Attempt temperature swaps on adjacent pairs with the given parity (odd/even sweep).
 */
static void pt_swap_pass(pt_sampler* s, int parity, double* inv_temps) {
    for (int i = parity; i < s->n_replicas - 1; i += 2) {
        double log_p = (inv_temps[i] - inv_temps[i + 1]) * (double)(s->h_energies[i] - s->h_energies[i + 1]);
        s->swap_attempts++;

        bool accept = (log_p >= 0.0) || (rng_uniform(&s->swap_rng) < exp(log_p));
        if (accept) {
            double tmp = s->h_temps[i];
            s->h_temps[i] = s->h_temps[i + 1];
            s->h_temps[i + 1] = tmp;

            double tmp_inv = inv_temps[i];
            inv_temps[i] = inv_temps[i + 1];
            inv_temps[i + 1] = tmp_inv;

            s->swap_accepted++;
        }
    }
}

/*
 * If any replica beat the best energy so far, copy that conformation into h_best.
 */
static void pt_update_best(pt_sampler* s) {
    for (int i = 0; i < s->n_replicas; i++) {
        if (s->h_energies[i] < s->h_best.energy) {
            // need to do ugly pointer arithmetic to get the ith chain on the host from GPU memory
            const char* src = reinterpret_cast<const char*>(s->d_chains) + (size_t)i * sizeof(hp_chain);
            CUDA_CHECK(cudaMemcpy(&s->h_best, src, sizeof(hp_chain), cudaMemcpyDeviceToHost));
        }
    }
}

extern "C" {

/*
 * Allocate replicas on the GPU, seed RNGs, and upload initial chains and temperatures.
 */
pt_sampler* pt_create(const char* seq, const double* temperatures, int n_replicas) {
    if (seq == NULL || temperatures == NULL || n_replicas < 2) return NULL;

    pt_sampler* s = (pt_sampler*)calloc(1, sizeof(pt_sampler));
    if (s == NULL) exit(EXIT_FAILURE);

    s->n_replicas = n_replicas;
    rng_seed(&s->swap_rng, 0xCAFEBABEL);

    size_t chain_bytes = (size_t)n_replicas * sizeof(hp_chain);
    hp_chain* host_chains = (hp_chain*)malloc(chain_bytes);
    if (host_chains == NULL) {
        pt_free(s);
        exit(EXIT_FAILURE);
    }

    s->h_temps = (double*)malloc((size_t)n_replicas * sizeof(double));
    s->h_energies = (int*)malloc((size_t)n_replicas * sizeof(int));
    if (s->h_temps == NULL || s->h_energies == NULL) {
        free(host_chains);
        pt_free(s);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < n_replicas; i++) {
        hp_chain* c = hp_chain_create(seq);
        if (c == NULL) {
            free(host_chains);
            pt_free(s);
            exit(EXIT_FAILURE);
        }
        host_chains[i] = *c;
        hp_chain_free(c);
        s->h_temps[i] = temperatures[i];
        s->h_energies[i] = host_chains[i].energy;
    }

    s->h_best = host_chains[0];
    for (int i = 1; i < n_replicas; i++) {
        if (host_chains[i].energy < s->h_best.energy) {
            s->h_best = host_chains[i];
        }
    }

    CUDA_CHECK(cudaMalloc(&s->d_chains, chain_bytes));
    CUDA_CHECK(cudaMalloc(&s->d_temps, (size_t)n_replicas * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&s->d_rngs, (size_t)n_replicas * sizeof(rng_state)));
    CUDA_CHECK(cudaMalloc(&s->d_energies, (size_t)n_replicas * sizeof(int)));

    CUDA_CHECK(cudaMemcpy(s->d_chains, host_chains, chain_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s->d_temps, s->h_temps, (size_t)n_replicas * sizeof(double), cudaMemcpyHostToDevice));

    rng_state* host_rngs = (rng_state*)malloc((size_t)n_replicas * sizeof(rng_state));
    if (host_rngs == NULL) {
        free(host_chains);
        pt_free(s);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n_replicas; i++) {
        rng_seed(&host_rngs[i], (uint64_t)(0x1234ULL + (uint64_t)i));
    }
    CUDA_CHECK(cudaMemcpy(s->d_rngs, host_rngs, (size_t)n_replicas * sizeof(rng_state), cudaMemcpyHostToDevice));

    free(host_chains);
    free(host_rngs);
    return s;
}

/*
 * Free GPU buffers, host arrays, and the sampler struct.
 */
void pt_free(pt_sampler* s) {
    if (s == NULL) return;
    if (s->d_chains) cudaFree(s->d_chains);
    if (s->d_temps) cudaFree(s->d_temps);
    if (s->d_rngs) cudaFree(s->d_rngs);
    if (s->d_energies) cudaFree(s->d_energies);
    free(s->h_energies);
    free(s->h_temps);
    free(s);
}

/*
 * Run n_swaps rounds: MH on the GPU, then host temperature swaps (odd/even), update best.
 */
void pt_run(pt_sampler* s, int steps_per_swap, int n_swaps) {
    if (s == NULL || steps_per_swap <= 0 || n_swaps <= 0) return;

    const int threads = 256;
    const int blocks = pt_launch_blocks(s->n_replicas);

    double* inv_temps = (double*)malloc((size_t)s->n_replicas * sizeof(double));
    if (inv_temps == NULL) exit(EXIT_FAILURE);

    for (int round = 0; round < n_swaps; round++) {
        pt_mh_kernel<<<blocks, threads>>>(s->d_chains, s->d_temps, s->d_rngs,
                                          s->n_replicas, steps_per_swap);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        pt_energies_kernel<<<blocks, threads>>>(s->d_chains, s->d_energies, s->n_replicas);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaMemcpy(s->h_energies, s->d_energies,
                             (size_t)s->n_replicas * sizeof(int), cudaMemcpyDeviceToHost));

        for (int i = 0; i < s->n_replicas; i++) {
            inv_temps[i] = 1.0 / s->h_temps[i];
        }

        pt_swap_pass(s, 0, inv_temps);
        pt_swap_pass(s, 1, inv_temps);

        CUDA_CHECK(cudaMemcpy(s->d_temps, s->h_temps,
                             (size_t)s->n_replicas * sizeof(double), cudaMemcpyHostToDevice));

        pt_update_best(s);
    }

    free(inv_temps);
}

/*
 * Pointer to the lowest-energy conformation seen during the run.
 */
const hp_chain* pt_best_chain(const pt_sampler* s) {
    if (s == NULL) return NULL;
    return &s->h_best;
}

/*
 * Fraction of attempted replica temperature swaps that were accepted.
 */
double pt_swap_acceptance_rate(const pt_sampler* s) {
    if (s == NULL || s->swap_attempts == 0) return 0.0;
    return (double)s->swap_accepted / (double)s->swap_attempts;
}

} /* extern "C" */
