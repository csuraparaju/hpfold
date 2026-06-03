#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "include/hp.h"
#include "include/pt.h"
#include "include/rng.h"

#define SOA_IDX(r, i) ((r) * HP_MAX_N + (i))

#define CUDA_CHECK(call)                                                       \
    do {                                                                       \
        cudaError_t err = (call);                                              \
        if (err != cudaSuccess) {                                              \
            fprintf(stderr, "CUDA error %s:%d: %s\n", __FILE__, __LINE__,      \
                    cudaGetErrorString(err));                                  \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
    } while (0)

__device__ int d_residue_at(const int* d_x, const int* d_y, int r, int n, int x, int y) {
    for (int i = 0; i < n; i++) {
        if (d_x[SOA_IDX(r, i)] == x && d_y[SOA_IDX(r, i)] == y) return i;
    }
    return -1;
}

__device__ uint64_t d_rng_u64(rng_state* rng) {
    uint64_t x = rng->state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    rng->state = x;
    return x;
}

__device__ double d_rng_uniform(rng_state* rng) {
    return (double)(d_rng_u64(rng) >> 11) * (1.0 / 9007199254740992.0);
}

__device__ int d_hp_chain_contacts_at(const int* d_x, const int* d_y, const char* d_aa,
                                      int r, int n, int aa_idx, int cx, int cy) {
    const int dx[4] = {1, -1, 0, 0};
    const int dy[4] = {0, 0, 1, -1};
    int contacts = 0;

    if (d_aa[aa_idx] != 'H') return 0;

    for (int d = 0; d < 4; d++) {
        int nx = cx + dx[d];
        int ny = cy + dy[d];
        if (nx < 0 || nx >= HP_LATTICE || ny < 0 || ny >= HP_LATTICE) continue;
        int j = d_residue_at(d_x, d_y, r, n, nx, ny);
        if (j == -1) continue;
        if (j == aa_idx) continue;
        if (j == aa_idx - 1 || j == aa_idx + 1) continue;
        if (d_aa[j] != 'H') continue;
        contacts++;
    }
    return contacts;
}

__device__ void d_hp_chain_get_coord(const int* d_x, const int* d_y, int r, int i,
                                     int* x, int* y) {
    *x = d_x[SOA_IDX(r, i)];
    *y = d_y[SOA_IDX(r, i)];
}

__device__ int d_hp_chain_site_occ(const int* d_x, const int* d_y, int r, int n, int x, int y) {
    return d_residue_at(d_x, d_y, r, n, x, y);
}

__device__ void d_hp_chain_commit_move(int* d_x, int* d_y, int* d_energy, const char* d_aa,
                                      int r, int n, int aa_idx,
                                      int old_x, int old_y, int new_x, int new_y) {
    int old_contacts = d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx, old_x, old_y);
    int new_contacts = d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx, new_x, new_y);
    int dE = -(new_contacts - old_contacts);

    d_x[SOA_IDX(r, aa_idx)] = new_x;
    d_y[SOA_IDX(r, aa_idx)] = new_y;
    d_energy[r] += dE;
}

__device__ void d_hp_chain_commit_two_site_move(int* d_x, int* d_y, int* d_energy, const char* d_aa,
    int r, int n,
    int aa_idx1, int old_x1, int old_y1, int new_x1, int new_y1,
    int aa_idx2, int old_x2, int old_y2, int new_x2, int new_y2) {
    int old_contacts = d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx1, old_x1, old_y1) +
                       d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx2, old_x2, old_y2);

    d_x[SOA_IDX(r, aa_idx1)] = new_x1;
    d_y[SOA_IDX(r, aa_idx1)] = new_y1;
    d_x[SOA_IDX(r, aa_idx2)] = new_x2;
    d_y[SOA_IDX(r, aa_idx2)] = new_y2;

    int new_contacts = d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx1, new_x1, new_y1) +
                       d_hp_chain_contacts_at(d_x, d_y, d_aa, r, n, aa_idx2, new_x2, new_y2);

    d_energy[r] += -(new_contacts - old_contacts);
}

__device__ bool d_in_bounds(int x, int y) {
    return x >= 0 && x < HP_LATTICE && y >= 0 && y < HP_LATTICE;
}

__device__ bool d_hp_end_flip(int* d_x, int* d_y, int* d_energy, const char* d_aa,
                              int r, int n, int aa_idx, rng_state* rng) {
    const int DX[4] = {1, -1, 0, 0};
    const int DY[4] = {0, 0, 1, -1};

    int nbor_idx = (aa_idx == 0) ? 1 : aa_idx - 1;
    int aa_x, aa_y, nbor_x, nbor_y;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx, &aa_x, &aa_y);
    d_hp_chain_get_coord(d_x, d_y, r, nbor_idx, &nbor_x, &nbor_y);

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
    if (d_hp_chain_site_occ(d_x, d_y, r, n, new_x, new_y) != -1) return false;

    d_hp_chain_commit_move(d_x, d_y, d_energy, d_aa, r, n, aa_idx, aa_x, aa_y, new_x, new_y);
    return true;
}

__device__ bool d_hp_corner_flip(int* d_x, int* d_y, int* d_energy, const char* d_aa,
                                 int r, int n, int aa_idx) {
    int prev_x, prev_y, next_x, next_y;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx - 1, &prev_x, &prev_y);
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx + 1, &next_x, &next_y);

    if (prev_x == next_x || prev_y == next_y) return false;

    int aa_x, aa_y;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx, &aa_x, &aa_y);

    int new_x = (aa_x == prev_x) ? next_x : prev_x;
    int new_y = (aa_y == prev_y) ? next_y : prev_y;

    if (new_x == aa_x && new_y == aa_y) return false;
    if (d_hp_chain_site_occ(d_x, d_y, r, n, new_x, new_y) != -1) return false;

    d_hp_chain_commit_move(d_x, d_y, d_energy, d_aa, r, n, aa_idx, aa_x, aa_y, new_x, new_y);
    return true;
}

__device__ bool d_hp_crankshaft(int* d_x, int* d_y, int* d_energy, const char* d_aa,
                                int r, int n, int aa_idx) {
    int ax, ay, bx, by;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx - 1, &ax, &ay);
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx + 2, &bx, &by);

    if (abs(ax - bx) + abs(ay - by) != 1) return false;

    int i_x, i_y, i1_x, i1_y;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx, &i_x, &i_y);
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx + 1, &i1_x, &i1_y);

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
    if (d_hp_chain_site_occ(d_x, d_y, r, n, new_i_x, new_i_y) != -1) return false;
    if (d_hp_chain_site_occ(d_x, d_y, r, n, new_i1_x, new_i1_y) != -1) return false;

    d_hp_chain_commit_two_site_move(d_x, d_y, d_energy, d_aa, r, n,
        aa_idx, i_x, i_y, new_i_x, new_i_y,
        aa_idx + 1, i1_x, i1_y, new_i1_x, new_i1_y);
    return true;
}

__device__ void d_mh_step(int* d_x, int* d_y, int* d_energy, const char* d_aa,
                          int r, int n, rng_state* rng, double temperature) {
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

    int old_energy = d_energy[r];
    int old_x1, old_y1, old_x2, old_y2;
    d_hp_chain_get_coord(d_x, d_y, r, aa_idx, &old_x1, &old_y1);
    if (move_type == 2) {
        d_hp_chain_get_coord(d_x, d_y, r, aa_idx + 1, &old_x2, &old_y2);
    }

    bool moved;
    if (move_type == 0) {
        moved = d_hp_end_flip(d_x, d_y, d_energy, d_aa, r, n, aa_idx, rng);
    } else if (move_type == 1) {
        moved = d_hp_corner_flip(d_x, d_y, d_energy, d_aa, r, n, aa_idx);
    } else {
        moved = d_hp_crankshaft(d_x, d_y, d_energy, d_aa, r, n, aa_idx);
    }

    if (!moved) return;

    int dE = d_energy[r] - old_energy;
    bool accept = (dE <= 0) || (d_rng_uniform(rng) < exp(-(double)dE / temperature));

    if (!accept) {
        int new_x1, new_y1;
        d_hp_chain_get_coord(d_x, d_y, r, aa_idx, &new_x1, &new_y1);
        if (move_type == 2) {
            int new_x2, new_y2;
            d_hp_chain_get_coord(d_x, d_y, r, aa_idx + 1, &new_x2, &new_y2);
            d_hp_chain_commit_two_site_move(d_x, d_y, d_energy, d_aa, r, n,
                aa_idx, new_x1, new_y1, old_x1, old_y1,
                aa_idx + 1, new_x2, new_y2, old_x2, old_y2);
        } else {
            d_hp_chain_commit_move(d_x, d_y, d_energy, d_aa, r, n,
                aa_idx, new_x1, new_y1, old_x1, old_y1);
        }
    }
}

/*
 * One thread per replica. State is SoA: d_x/d_y are [replica][residue], d_aa is shared.
 */
__global__ void pt_mh_kernel(int* d_x, int* d_y, const char* d_aa, int* d_energy,
                             double* d_temps, rng_state* d_rngs,
                             int chain_n, int n_replicas, int steps) {
    int r = blockIdx.x * blockDim.x + threadIdx.x;
    if (r >= n_replicas) return;

    double temp = d_temps[r];
    for (int s = 0; s < steps; s++) {
        d_mh_step(d_x, d_y, d_energy, d_aa, r, chain_n, &d_rngs[r], temp);
    }
}

struct pt_sampler {
    int* d_x;
    int* d_y;
    char* d_aa;
    int* d_energy;
    double* d_temps;
    rng_state* d_rngs;

    int* h_energies;
    double* h_temps;

    hp_chain* h_best;
    int* gather_x;
    int* gather_y;
    int chain_n;
    int n_replicas;
    long swap_attempts;
    long swap_accepted;

    rng_state swap_rng;
};

static int pt_launch_blocks(int n) {
    const int threads = 256;
    return (n + threads - 1) / threads;
}

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

static void pt_gather_coords(const pt_sampler* s, int replica, int* xs, int* ys) {
    size_t coord_bytes = (size_t)s->chain_n * sizeof(int);
    CUDA_CHECK(cudaMemcpy(xs, s->d_x + SOA_IDX(replica, 0), coord_bytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(ys, s->d_y + SOA_IDX(replica, 0), coord_bytes, cudaMemcpyDeviceToHost));
}

static void pt_update_best(pt_sampler* s) {
    for (int i = 0; i < s->n_replicas; i++) {
        if (s->h_energies[i] < hp_chain_energy(s->h_best)) {
            pt_gather_coords(s, i, s->gather_x, s->gather_y);
            hp_chain_set_coords(s->h_best, s->gather_x, s->gather_y, s->chain_n);
        }
    }
}

extern "C" {

pt_sampler* pt_create(const char* seq, const double* temperatures, int n_replicas) {
    if (seq == NULL || temperatures == NULL || n_replicas < 2) return NULL;

    pt_sampler* s = (pt_sampler*)calloc(1, sizeof(pt_sampler));
    if (s == NULL) exit(EXIT_FAILURE);

    s->n_replicas = n_replicas;
    rng_seed(&s->swap_rng, 0xCAFEBABEL);

    hp_chain* init = hp_chain_create(seq);
    s->chain_n = hp_chain_length(init);

    s->h_best = hp_chain_clone(init);

    size_t replica_slots = (size_t)n_replicas * HP_MAX_N;
    int* host_x = (int*)malloc(replica_slots * sizeof(int));
    int* host_y = (int*)malloc(replica_slots * sizeof(int));
    int* host_energy = (int*)malloc((size_t)n_replicas * sizeof(int));
    char* host_aa = (char*)malloc((size_t)s->chain_n);
    int* init_x = (int*)malloc((size_t)s->chain_n * sizeof(int));
    int* init_y = (int*)malloc((size_t)s->chain_n * sizeof(int));
    if (host_x == NULL || host_y == NULL || host_energy == NULL || host_aa == NULL ||
        init_x == NULL || init_y == NULL) {
        hp_chain_free(init);
        free(host_x);
        free(host_y);
        free(host_energy);
        free(host_aa);
        free(init_x);
        free(init_y);
        pt_free(s);
        exit(EXIT_FAILURE);
    }

    s->h_temps = (double*)malloc((size_t)n_replicas * sizeof(double));
    s->h_energies = (int*)malloc((size_t)n_replicas * sizeof(int));
    s->gather_x = (int*)malloc((size_t)s->chain_n * sizeof(int));
    s->gather_y = (int*)malloc((size_t)s->chain_n * sizeof(int));
    if (s->h_temps == NULL || s->h_energies == NULL || s->gather_x == NULL || s->gather_y == NULL) {
        hp_chain_free(init);
        free(host_x);
        free(host_y);
        free(host_energy);
        free(host_aa);
        free(init_x);
        free(init_y);
        pt_free(s);
        exit(EXIT_FAILURE);
    }

    hp_chain_get_coords(init, init_x, init_y, s->chain_n);
    int init_energy = hp_chain_energy(init);
    for (int j = 0; j < s->chain_n; j++) {
        host_aa[j] = hp_chain_aa_type(init, j);
    }

    for (int i = 0; i < n_replicas; i++) {
        for (int j = 0; j < s->chain_n; j++) {
            host_x[SOA_IDX(i, j)] = init_x[j];
            host_y[SOA_IDX(i, j)] = init_y[j];
        }
        host_energy[i] = init_energy;
        s->h_temps[i] = temperatures[i];
        s->h_energies[i] = init_energy;
    }

    CUDA_CHECK(cudaMalloc(&s->d_x, replica_slots * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&s->d_y, replica_slots * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&s->d_aa, HP_MAX_N * sizeof(char)));
    CUDA_CHECK(cudaMalloc(&s->d_energy, (size_t)n_replicas * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&s->d_temps, (size_t)n_replicas * sizeof(double)));
    CUDA_CHECK(cudaMalloc(&s->d_rngs, (size_t)n_replicas * sizeof(rng_state)));

    CUDA_CHECK(cudaMemcpy(s->d_x, host_x, replica_slots * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s->d_y, host_y, replica_slots * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemset(s->d_aa, 0, HP_MAX_N * sizeof(char)));
    CUDA_CHECK(cudaMemcpy(s->d_aa, host_aa, (size_t)s->chain_n, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s->d_energy, host_energy, (size_t)n_replicas * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(s->d_temps, s->h_temps, (size_t)n_replicas * sizeof(double), cudaMemcpyHostToDevice));

    rng_state* host_rngs = (rng_state*)malloc((size_t)n_replicas * sizeof(rng_state));
    if (host_rngs == NULL) {
        hp_chain_free(init);
        free(host_x);
        free(host_y);
        free(host_energy);
        free(host_aa);
        free(init_x);
        free(init_y);
        pt_free(s);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n_replicas; i++) {
        rng_seed(&host_rngs[i], (uint64_t)(0x1234ULL + (uint64_t)i));
    }
    CUDA_CHECK(cudaMemcpy(s->d_rngs, host_rngs, (size_t)n_replicas * sizeof(rng_state), cudaMemcpyHostToDevice));

    hp_chain_free(init);
    free(host_x);
    free(host_y);
    free(host_energy);
    free(host_aa);
    free(init_x);
    free(init_y);
    free(host_rngs);
    return s;
}

void pt_free(pt_sampler* s) {
    if (s == NULL) return;
    if (s->d_x) cudaFree(s->d_x);
    if (s->d_y) cudaFree(s->d_y);
    if (s->d_aa) cudaFree(s->d_aa);
    if (s->d_energy) cudaFree(s->d_energy);
    if (s->d_temps) cudaFree(s->d_temps);
    if (s->d_rngs) cudaFree(s->d_rngs);
    if (s->h_best) hp_chain_free(s->h_best);
    free(s->gather_x);
    free(s->gather_y);
    free(s->h_energies);
    free(s->h_temps);
    free(s);
}

void pt_run(pt_sampler* s, int steps_per_swap, int n_swaps) {
    if (s == NULL || steps_per_swap <= 0 || n_swaps <= 0) return;

    const int threads = 256;
    const int blocks = pt_launch_blocks(s->n_replicas);

    double* inv_temps = (double*)malloc((size_t)s->n_replicas * sizeof(double));
    if (inv_temps == NULL) exit(EXIT_FAILURE);

    for (int round = 0; round < n_swaps; round++) {
        pt_mh_kernel<<<blocks, threads>>>(s->d_x, s->d_y, s->d_aa, s->d_energy,
                                          s->d_temps, s->d_rngs,
                                          s->chain_n, s->n_replicas, steps_per_swap);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaDeviceSynchronize());

        CUDA_CHECK(cudaMemcpy(s->h_energies, s->d_energy,
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

const hp_chain* pt_best_chain(const pt_sampler* s) {
    if (s == NULL) return NULL;
    return s->h_best;
}

double pt_swap_acceptance_rate(const pt_sampler* s) {
    if (s == NULL || s->swap_attempts == 0) return 0.0;
    return (double)s->swap_accepted / (double)s->swap_attempts;
}

} /* extern "C" */
