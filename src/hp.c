#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "include/contracts.h"
#include "include/hp.h"

#define ORIGIN      (HP_LATTICE / 2)
#define LATTICE_IDX(x, y) ((x) * HP_LATTICE + (y))

struct hp_chain {
    int x[HP_MAX_N];
    int y[HP_MAX_N];
    int n;
    char aa_seq[HP_MAX_N];
    int lattice_occ[HP_LATTICE * HP_LATTICE];
    int energy;
};

int hp_chain_contacts_at(const hp_chain* chain, int aa_idx, int cx, int cy) {
    REQUIRES(chain != NULL);
    REQUIRES(aa_idx >= 0 && aa_idx < chain->n);
    REQUIRES(cx >= 0 && cx < HP_LATTICE && cy >= 0 && cy < HP_LATTICE);

    static const int dx[4] = {1, -1, 0, 0};
    static const int dy[4] = {0, 0, 1, -1};
    int contacts = 0;

    if (chain->aa_seq[aa_idx] != 'H') return 0;

    for (int d = 0; d < 4; d++) {
        int nx = cx + dx[d];
        int ny = cy + dy[d];
        if (nx < 0 || nx >= HP_LATTICE || ny < 0 || ny >= HP_LATTICE) continue;
        int j = chain->lattice_occ[LATTICE_IDX(nx, ny)];

        if (j == -1) continue;
        if (j == aa_idx) continue;
        if (j == aa_idx - 1 || j == aa_idx + 1) continue;
        if (chain->aa_seq[j] != 'H') continue;

        contacts++;
    }
    return contacts;
}

static int compute_energy(const hp_chain* chain) {
    REQUIRES(chain != NULL);
    REQUIRES(chain->n > 0 && chain->n <= HP_MAX_N);
    int contacts = 0;

    for (int i = 0; i < chain->n; i++) {
        contacts += hp_chain_contacts_at(chain, i, chain->x[i], chain->y[i]);
    }
    return -contacts / 2;
}

bool is_hp_chain(const hp_chain* chain) {
    if (chain == NULL) return false;
    if (chain->n <= 0 || chain->n > HP_MAX_N) return false;

    int num_occupied = 0;
    for (int idx = 0; idx < HP_LATTICE * HP_LATTICE; idx++) {
        int b = chain->lattice_occ[idx];
        if (b == -1) continue;
        if (b < 0 || b >= chain->n) return false;
        num_occupied++;
    }
    if (num_occupied != chain->n) return false;

    for (int i = 0; i < chain->n; i++) {
        if (chain->x[i] < 0 || chain->x[i] >= HP_LATTICE) return false;
        if (chain->y[i] < 0 || chain->y[i] >= HP_LATTICE) return false;
        if (chain->aa_seq[i] != 'H' && chain->aa_seq[i] != 'P') return false;

        int site = LATTICE_IDX(chain->x[i], chain->y[i]);
        if (chain->lattice_occ[site] != i) return false;
    }

    for (int i = 0; i < chain->n - 1; i++) {
        int mdist = abs(chain->x[i] - chain->x[i+1]) + abs(chain->y[i] - chain->y[i+1]);
        if (mdist != 1) return false;
    }

    if (chain->energy != compute_energy(chain)) return false;

    return true;
}

hp_chain* hp_chain_create(const char* aa_seq) {
    REQUIRES(aa_seq != NULL);
    int n = (int)strlen(aa_seq);
    REQUIRES(n > 0 && n <= HP_MAX_N);

    hp_chain* chain = malloc(sizeof(hp_chain));
    if (chain == NULL) exit(EXIT_FAILURE);

    chain->n = n;

    for (int i = 0; i < HP_LATTICE * HP_LATTICE; i++) chain->lattice_occ[i] = -1;

    for (int i = 0; i < n; i++) {
        chain->aa_seq[i] = aa_seq[i];
        chain->x[i] = ORIGIN;
        chain->y[i] = ORIGIN + i;
        chain->lattice_occ[LATTICE_IDX(ORIGIN, ORIGIN + i)] = i;
    }

    chain->energy = compute_energy(chain);
    ASSERT(chain->energy == 0);

    ENSURES(is_hp_chain(chain));
    return chain;
}

void hp_chain_free(hp_chain* chain) {
    REQUIRES(is_hp_chain(chain));
    free(chain);
}

int hp_chain_length(const hp_chain* chain) {
    REQUIRES(is_hp_chain(chain));
    return chain->n;
}

char hp_chain_aa_type(const hp_chain* chain, int i) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(i >= 0 && i < chain->n);
    return chain->aa_seq[i];
}

int hp_chain_energy(const hp_chain* chain) {
    REQUIRES(is_hp_chain(chain));
    return chain->energy;
}

void hp_chain_get_coords(const hp_chain* chain, int* xs, int* ys, int n) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(xs != NULL && ys != NULL);
    REQUIRES(n == chain->n);
    for (int i = 0; i < n; i++) {
        xs[i] = chain->x[i];
        ys[i] = chain->y[i];
    }
}

void hp_chain_get_coord(const hp_chain* chain, int i, int* x, int* y) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(i >= 0 && i < chain->n);
    REQUIRES(x != NULL && y != NULL);
    *x = chain->x[i];
    *y = chain->y[i];
}

int hp_chain_site_occ(const hp_chain* chain, int x, int y) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(x >= 0 && x < HP_LATTICE && y >= 0 && y < HP_LATTICE);
    return chain->lattice_occ[LATTICE_IDX(x, y)];
}

void hp_chain_commit_move(hp_chain* chain,
    int aa_idx,
    int old_x, int old_y,
    int new_x, int new_y
) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(aa_idx >= 0 && aa_idx < chain->n);
    REQUIRES(old_x >= 0 && old_x < HP_LATTICE && old_y >= 0 && old_y < HP_LATTICE);
    REQUIRES(new_x >= 0 && new_x < HP_LATTICE && new_y >= 0 && new_y < HP_LATTICE);

    int old_contacts = hp_chain_contacts_at(chain, aa_idx, old_x, old_y);
    int new_contacts = hp_chain_contacts_at(chain, aa_idx, new_x, new_y);
    int dE = -(new_contacts - old_contacts);

    chain->lattice_occ[LATTICE_IDX(old_x, old_y)] = -1;
    chain->lattice_occ[LATTICE_IDX(new_x, new_y)] = aa_idx;
    chain->x[aa_idx] = new_x;
    chain->y[aa_idx] = new_y;
    chain->energy += dE;
}

void hp_chain_commit_two_site_move(hp_chain* chain,
    int aa_idx1, int old_x1, int old_y1, int new_x1, int new_y1,
    int aa_idx2, int old_x2, int old_y2, int new_x2, int new_y2
) {
    REQUIRES(is_hp_chain(chain));
    REQUIRES(aa_idx1 >= 0 && aa_idx1 < chain->n);
    REQUIRES(aa_idx2 >= 0 && aa_idx2 < chain->n);
    REQUIRES(old_x1 >= 0 && old_x1 < HP_LATTICE && old_y1 >= 0 && old_y1 < HP_LATTICE);
    REQUIRES(old_x2 >= 0 && old_x2 < HP_LATTICE && old_y2 >= 0 && old_y2 < HP_LATTICE);
    REQUIRES(new_x1 >= 0 && new_x1 < HP_LATTICE && new_y1 >= 0 && new_y1 < HP_LATTICE);
    REQUIRES(new_x2 >= 0 && new_x2 < HP_LATTICE && new_y2 >= 0 && new_y2 < HP_LATTICE);

    // Compute old contacts before clearing any sites
    int old_contacts = hp_chain_contacts_at(chain, aa_idx1, old_x1, old_y1) +
                       hp_chain_contacts_at(chain, aa_idx2, old_x2, old_y2);

    // Clear old positions then place at new positions
    chain->lattice_occ[LATTICE_IDX(old_x1, old_y1)] = -1;
    chain->lattice_occ[LATTICE_IDX(old_x2, old_y2)] = -1;
    chain->lattice_occ[LATTICE_IDX(new_x1, new_y1)] = aa_idx1;
    chain->lattice_occ[LATTICE_IDX(new_x2, new_y2)] = aa_idx2;

    chain->x[aa_idx1] = new_x1;
    chain->y[aa_idx1] = new_y1;
    chain->x[aa_idx2] = new_x2;
    chain->y[aa_idx2] = new_y2;

    // Compute new contacts after the move so neighbours reflect the new state
    int new_contacts = hp_chain_contacts_at(chain, aa_idx1, new_x1, new_y1) +
                       hp_chain_contacts_at(chain, aa_idx2, new_x2, new_y2);

    chain->energy += -(new_contacts - old_contacts);
}

void hp_chain_print(const hp_chain* chain) {
    REQUIRES(is_hp_chain(chain));

    const int ROWS = 2 * HP_LATTICE - 1;
    const int COLS = 4 * HP_LATTICE - 1;

    char canvas[ROWS][COLS + 1];

    for (int r = 0; r < ROWS; r++) {
        for (int c = 0; c < COLS; c++) canvas[r][c] = ' ';
        canvas[r][COLS] = '\0';
    }

    for (int i = 0; i < chain->n; i++) {
        int x = chain->x[i];
        int y = chain->y[i];

        int row = 2 * x;
        int col = 4 * y;

        char label[8];
        snprintf(label, sizeof(label), "%c", chain->aa_seq[i]);

        for (int k = 0; label[k] != '\0' && col + k < COLS; k++)
            canvas[row][col + k] = label[k];

        if (i < chain->n - 1) {
            int nx = chain->x[i + 1];
            int ny = chain->y[i + 1];

            int drow = 2 * (nx - x);
            int dcol = 4 * (ny - y);

            if (drow != 0) {
                canvas[row + drow / 2][col] = '|';
            }

            if (dcol != 0) {
                int step = (dcol > 0) ? 1 : -1;
                for (int c = col + step; c != col + dcol; c += step)
                    canvas[row][c] = '-';
            }
        }
    }

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

    printf("\nEnergy = %d\n", chain->energy);

    for (int r = min_r; r <= max_r; r++) {
        for (int c = min_c; c <= max_c; c++) putchar(canvas[r][c]);
        putchar('\n');
    }
    putchar('\n');
}

hp_chain* hp_chain_clone(const hp_chain* chain) {
    REQUIRES(is_hp_chain(chain));
    hp_chain* out = malloc(sizeof(hp_chain));
    *out = *chain;
    ENSURES(is_hp_chain(out));
    return out;
}
