#pragma once
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif


#define HP_MAX_N   64
#define HP_LATTICE (2 * HP_MAX_N)

typedef struct hp_chain hp_chain;

/**
 * Allocate and initialize an HP chain representation for
 * the given null-terminated HP sequence (e.g. "HPPHH").
 * This representation contains information about the folded
 * structure of the chain, as given on the 2D lattice.
 */
hp_chain* hp_chain_create(const char* seq);

/**
 * Free all memory associated with this HP chain.
 */
void hp_chain_free(hp_chain* c);


/**
 * Return the number of amino acids in the chain.
 */
int hp_chain_length(const hp_chain* c);

/**
 * Return the amino acid type ('H' or 'P') at position i.
 */
char hp_chain_aa_type(const hp_chain* c, int i);

/**
 * Return the energy associated with this chain.
 */
int hp_chain_energy(const hp_chain* c);

/**
 * Get the 2D lattice coordinates of each amino acid in the chain's folded
 * configuration. Writes into xs and ys such that xs[i], ys[i] is the
 * position of the ith amino acid. n must equal hp_chain_length(c).
 */
void hp_chain_get_coords(const hp_chain* c, int* xs, int* ys, int n);

/**
 * Get the 2D lattice coordinates of a single amino acid i.
 */
void hp_chain_get_coord(const hp_chain* c, int i, int* x, int* y);

/**
 * Return the amino acid index occupying lattice site (x, y), or -1 if empty.
 */
int hp_chain_site_occ(const hp_chain* c, int x, int y);

/**
 * Return the number of non-bonded H-H contacts that amino acid aa_idx would
 * have if placed at (cx, cy). Can be used to compute energy changes incrementally.
 */
int hp_chain_contacts_at(const hp_chain* c, int aa_idx, int cx, int cy);

/**
 * Move amino acid aa_idx from (old_x, old_y) to (new_x, new_y), updating
 * coordinates and energy. The caller must have verified that the destination
 * is in-bounds and empty.
 */
void hp_chain_commit_move(hp_chain* c, int aa_idx, int old_x, int old_y, int new_x, int new_y);

/**
 * Move two amino acids atomically, computing the combined energy delta and
 * applying all updates in the correct order. Can be used for crankshaft moves.
 * The caller must have verified that all destination sites are in-bounds and empty.
 */
void hp_chain_commit_two_site_move(hp_chain* c,
                                   int aa_idx1, int old_x1, int old_y1,
                                   int new_x1,  int new_y1,
                                   int aa_idx2, int old_x2, int old_y2,
                                   int new_x2,  int new_y2);

/**
 * Pretty print the chain.
 */
void hp_chain_print(const hp_chain* c);

/**
 * Deep copy of the chain.
 */
hp_chain* hp_chain_clone(const hp_chain* c);

#ifdef __cplusplus
}
#endif
