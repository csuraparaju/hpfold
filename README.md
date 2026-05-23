# hpfold

Monte Carlo folding of HP lattice proteins on a 2D square lattice, using Metropolis–Hastings sampling with local conformation moves.

## Overview

The [HP model](https://en.wikipedia.org/wiki/Hydrophobic-polar_protein_folding_model) simplifies a protein as a chain of H (hydrophobic) and P (polar) residues on a 2D lattice. The chain must stay self-avoiding with consecutive residues on adjacent lattice sites. Energy comes only from non-bonded H–H contacts (neighbors in 2D space that are not neighbors along the sequence).

This project provides:

- A `hp_chain` representation that tracks residue coordinates on the 2D lattice with incremental energy updates (occupancy by coordinate scan, no occupancy grid)
- Three local move types (end flip, corner flip, crankshaft)
- A Metropolis–Hastings sampler over conformations at a fixed temperature
- Optional **CUDA parallel tempering** (`mc_hp_pt`): one GPU thread per replica, temperature swaps on the host

## Moves

| Move | Description |
|------|-------------|
| **End flip** | An endpoint pivots around its sole sequence neighbor; one of three lattice neighbors is chosen at random. |
| **Corner flip** | An interior residue at a "corner" of a unit square formed by its sequence neighbors is reflected to the opposite corner. |
| **Crankshaft** | Two consecutive interior residues rotate when the residues two steps away form a U-shaped anchor on the lattice. |

The MH sampler picks a move type uniformly, proposes a valid residue index for that type, applies the move, then accepts or rejects using the standard criterion with temperature `T`:

`accept if dE <= 0`, else accept with probability `exp(-dE / T)`.

Rejected moves are undone by committing the reverse displacement.

## Project layout

    hpfold/
    ├── main.c              # Demo: MH run on a fixed HP sequence
    ├── makefile            # Build mc_hp (release) and mc_hp_dbg (debug)
    ├── include/
    │   ├── hp.h            # HP chain API (coords, energy, moves)
    │   ├── mh.h            # Metropolis–Hastings sampler API
    │   ├── rng.h           # xorshift64 RNG
    │   ├── contracts.h     # REQUIRES/ENSURES/ASSERT (15-122 style)
    │   └── pt.h            # Parallel tempering sampler API (CUDA)
    └── src/
        ├── hp.c            # Chain representation and energy
        ├── moves.c         # End flip, corner flip, crankshaft
        ├── mh.c            # MH step/run and best-chain tracking
        ├── rng.c           # RNG implementation
        └── pt.cu           # CUDA parallel tempering (device MH + host swaps)

## Requirements

- **GCC** (or compatible C compiler)
- **C99**
- **libm** (`-lm`)
- **CUDA** (`nvcc`) for `mc_hp_pt` only

## Build and run

```bash
make          # builds ./mc_hp
make mc_hp_pt # builds ./mc_hp_pt (requires nvcc + CUDA GPU)
make debug    # builds ./mc_hp_dbg with -g -DDEBUG (contracts enabled)
make clean    # removes binaries
```

Run the demo (defaults match the original hard-coded values):

```bash
./mc_hp
./mc_hp -h
./mc_hp -s HPPHPPHHPPHHPPHH -t 5.0 -n 50000
```

`main.c` options:

| Flag | Meaning | Default |
|------|---------|---------|
| `-s` | HP sequence (`H` / `P`) | `HPPHPPHHPPHHPPHH` |
| `-t` | Temperature | `5.0` |
| `-n` | MH steps | `50000` |

The program prints:

- Initial and final conformations (ASCII lattice drawing)
- Best conformation seen during the run
- Best energy and overall acceptance rate

### Parallel tempering (CUDA)

```bash
make mc_hp_pt
./mc_hp_pt
./mc_hp_pt -h
./mc_hp_pt -s HPPHPPHHPPHHPPHH -r 8 -l 1 -u 10 -m 5000 -n 200
```

`main_pt.c` options:

| Flag | Meaning | Default |
|------|---------|---------|
| `-s` | HP sequence | `HPPHPPHHPPHHPPHH` |
| `-r` | Number of replicas | `8` |
| `-l` | Minimum temperature | `1.0` |
| `-u` | Maximum temperature | `10.0` |
| `-m` | MH steps per replica between swaps | `5000` |
| `-n` | Swap rounds | `200` |

Replicas use a geometric temperature ladder from `-l` to `-u`. Swaps exchange **temperatures** (not 65 KB chain copies); only replica energies are copied back for Metropolis swap tests.

