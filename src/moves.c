#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "include/contracts.h"
#include "include/hp.h"
#include "include/rng.h"

/* Client code of hp interface for performing moves on an
 * HP chain. Implements three types of moves: end flip move,
 * corner flip move, and a crank shaft move
*/

static const int DX[4] = { 1, -1,  0,  0 };
static const int DY[4] = { 0,  0,  1, -1 };

static bool in_bounds(int x, int y) {
    return (x >= 0 && x < HP_LATTICE && y >= 0 && y < HP_LATTICE);
}

/*
 * Attempt an end-flip move on amino acid aa_idx, which must be one of the
 * two chain endpoints (0 or n-1).
 *
 * For an endpoint, there are exactly three candidate positions around the
 * adjacent residue. One is the current endpoint position, and the other two
 * are possible flipped positions.
 *
 * This function randomly chooses one of those three candidates, and commits the
 * move and updates energy iff it is valid.
 *
 * Returns true if the move was applied, false otherwise.
 */
bool hp_end_flip(hp_chain* chain, int aa_idx, rng_state* rng) {
    REQUIRES(chain != NULL);
    REQUIRES(rng != NULL);
    REQUIRES(aa_idx == 0 || aa_idx == hp_chain_length(chain) - 1);

    // Identify the adjacent residue the endpoint pivots around
    int nbor_idx = (aa_idx == 0) ? 1 : aa_idx - 1;

    int aa_x, aa_y, nbor_x, nbor_y;
    hp_chain_get_coord(chain, aa_idx,   &aa_x,   &aa_y);
    hp_chain_get_coord(chain, nbor_idx, &nbor_x, &nbor_y);

    ASSERT(abs(aa_x - nbor_x) + abs(aa_y - nbor_y) == 1);


    // Enumerate the four lattice neighbours of the adjacent residue, excluding
    // the current endpoint position
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

    ASSERT(n_cand == 3);

    int choice = (int)(rng_uniform(rng) * 3.0);
    if (choice < 0) choice = 0;
    if (choice > 2) choice = 2;

    int new_x = cand_x[choice];
    int new_y = cand_y[choice];

    if (!in_bounds(new_x, new_y)) return false;
    if (hp_chain_site_occ(chain, new_x, new_y) != -1) return false;

    hp_chain_commit_move(chain, aa_idx, aa_x, aa_y, new_x, new_y);

    return true;
}

/*
 * Attempt a corner-flip move on an interior amino acid aa_idx.
 *
 * This move is eligible when the residue sits at a corner of a unit square
 * formed by its two sequence neighbors. The residue is reflected to the
 * opposite corner of that square.
 *
 * Returns true if the move was applied, false otherwise.
 */
bool hp_corner_flip(hp_chain* chain, int aa_idx) {
    REQUIRES(chain != NULL);
    REQUIRES(aa_idx > 0 && aa_idx < hp_chain_length(chain) - 1);

    int prev_x, prev_y, next_x, next_y;
    hp_chain_get_coord(chain, aa_idx - 1, &prev_x, &prev_y);
    hp_chain_get_coord(chain, aa_idx + 1, &next_x, &next_y);

    // Corner flip is only possible if the two sequence neighbors are not
    // aligned in the same row or column
    if (prev_x == next_x || prev_y == next_y) return false;

    int aa_x, aa_y;
    hp_chain_get_coord(chain, aa_idx, &aa_x, &aa_y);


    // The residue sits at one corner of the square; the opposite corner is
    // determined by swapping the shared coordinate with each neighbor
    int new_x = (aa_x == prev_x) ? next_x : prev_x;
    int new_y = (aa_y == prev_y) ? next_y : prev_y;

    if (new_x == aa_x && new_y == aa_y) return false;

    ASSERT(in_bounds(new_x, new_y));

    if (hp_chain_site_occ(chain, new_x, new_y) != -1) return false;

    hp_chain_commit_move(chain, aa_idx, aa_x, aa_y, new_x, new_y);

    return true;
}

/*
 * Attempt a crankshaft move on the pair (aa_idx, aa_idx + 1).
 *
 * This move is only applied if the following conditions hold
 * 1. aa_idx must satisfy 1 <= aa_idx <= n-3
 * 2. the "anchor" residues aa_idx-1 and aa_idx+2 must be lattice-adjacent,
 *     forming a U-shape that aa_idx and aa_idx+1 hang off of.
 * 3. the two middle residues are reflected across the line connecting the anchors.
 *
 * Returns true if the move was applied, false otherwise.
 */
bool hp_crankshaft(hp_chain* chain, int aa_idx) {
    REQUIRES(chain != NULL);
    REQUIRES(aa_idx >= 1 && aa_idx <= hp_chain_length(chain) - 3);

    int ax, ay, bx, by;
    hp_chain_get_coord(chain, aa_idx - 1, &ax, &ay);
    hp_chain_get_coord(chain, aa_idx + 2, &bx, &by);


    // The two anchor residues must be lattice-adjacent. If not, the local
    // geometry does not form a crankshaft U-shape
    if (abs(ax - bx) + abs(ay - by) != 1) return false;

    int i_x, i_y, i1_x, i1_y;
    hp_chain_get_coord(chain, aa_idx,     &i_x,  &i_y);
    hp_chain_get_coord(chain, aa_idx + 1, &i1_x, &i1_y);


    // Reflect the two moving residues across the line connecting the anchors.
    // Since the anchors are adjacent, they share either their x or y coordinate
    int new_i_x, new_i_y, new_i1_x, new_i1_y;

    if (ax == bx) {
        new_i_x  = 2 * ax - i_x;
        new_i_y  = i_y;

        new_i1_x = 2 * ax - i1_x;
        new_i1_y = i1_y;
    } else {
        ASSERT(ay == by);
        new_i_x  = i_x;
        new_i_y  = 2 * ay - i_y;

        new_i1_x = i1_x;
        new_i1_y = 2 * ay - i1_y;
    }

    if (new_i_x == i_x && new_i_y == i_y &&
        new_i1_x == i1_x && new_i1_y == i1_y) return false;

    if (!in_bounds(new_i_x, new_i_y) || !in_bounds(new_i1_x, new_i1_y)) return false;

    if (hp_chain_site_occ(chain, new_i_x,  new_i_y)  != -1) return false;
    if (hp_chain_site_occ(chain, new_i1_x, new_i1_y) != -1) return false;

    hp_chain_commit_two_site_move(chain,
        aa_idx,     i_x,  i_y,  new_i_x,  new_i_y,
        aa_idx + 1, i1_x, i1_y, new_i1_x, new_i1_y);

    return true;
}
