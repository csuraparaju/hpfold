#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "include/hp.h"
#include "include/pt.h"

static int valid_hp_sequence(const char* seq) {
    if (seq == NULL || seq[0] == '\0') {
        return 0;
    }
    if (strlen(seq) > (size_t)HP_MAX_N) {
        return 0;
    }
    for (size_t i = 0; seq[i] != '\0'; i++) {
        if (seq[i] != 'H' && seq[i] != 'P') {
            return 0;
        }
    }
    return 1;
}

static void fill_temperature_ladder(double* temps, int n, double t_min, double t_max) {
    if (n == 1) {
        temps[0] = t_min;
        return;
    }
    for (int i = 0; i < n; i++) {
        double frac = (double)i / (double)(n - 1);
        temps[i] = t_min * pow(t_max / t_min, frac);
    }
}

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [-s sequence] [-r replicas] [-l t_min] [-u t_max]\n"
        "              [-m steps_per_swap] [-n swaps]\n"
        "\n"
        "  -s sequence        HP amino acid string (H/P), max length %d\n"
        "  -r replicas        Number of parallel replicas (>= 2), default 8\n"
        "  -l t_min           Lowest temperature (> 0), default 1.0\n"
        "  -u t_max           Highest temperature (> t_min), default 10.0\n"
        "  -m steps_per_swap  MH steps per replica between swaps, default 5000\n"
        "  -n swaps           Number of swap rounds, default 200\n"
        "  -h                 Show this help\n"
        "\n"
        "Temperatures form a geometric ladder from t_min to t_max.\n"
        "\n"
        "Example:\n"
        "  %s -s HPPHPPHHPPHHPPHH -r 8 -l 1 -u 10 -m 5000 -n 200\n",
        prog, HP_MAX_N, prog);
}

int main(int argc, char** argv) {
    const char* seq = "HPPHPPHHPPHHPPHH";
    int n_replicas = 8;
    double t_min = 1.0;
    double t_max = 10.0;
    int steps_per_swap = 5000;
    int n_swaps = 200;

    int opt;
    while ((opt = getopt(argc, argv, "s:r:l:u:m:n:h")) != -1) {
        switch (opt) {
            case 's':
                seq = optarg;
                break;
            case 'r':
                n_replicas = atoi(optarg);
                break;
            case 'l':
                t_min = atof(optarg);
                break;
            case 'u':
                t_max = atof(optarg);
                break;
            case 'm':
                steps_per_swap = atoi(optarg);
                break;
            case 'n':
                n_swaps = atoi(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (!valid_hp_sequence(seq)) {
        fprintf(stderr, "error: invalid sequence (use H/P, length 1..%d)\n", HP_MAX_N);
        return 1;
    }
    if (n_replicas < 2) {
        fprintf(stderr, "error: replicas must be >= 2\n");
        return 1;
    }
    if (t_min <= 0.0 || t_max <= t_min) {
        fprintf(stderr, "error: require 0 < t_min < t_max\n");
        return 1;
    }
    if (steps_per_swap <= 0 || n_swaps <= 0) {
        fprintf(stderr, "error: steps_per_swap and swaps must be > 0\n");
        return 1;
    }

    double* temperatures = (double*)malloc((size_t)n_replicas * sizeof(double));
    if (temperatures == NULL) {
        fprintf(stderr, "error: out of memory\n");
        return 1;
    }
    fill_temperature_ladder(temperatures, n_replicas, t_min, t_max);

    pt_sampler* s = pt_create(seq, temperatures, n_replicas);
    free(temperatures);

    if (s == NULL) {
        fprintf(stderr, "error: pt_create failed\n");
        return 1;
    }

    printf("sequence:         %s\n", seq);
    printf("replicas:         %d\n", n_replicas);
    printf("temperatures:     %.3f .. %.3f (geometric)\n", t_min, t_max);
    printf("steps_per_swap:   %d\n", steps_per_swap);
    printf("swap rounds:      %d\n\n", n_swaps);

    pt_run(s, steps_per_swap, n_swaps);

    printf("Best chain (energy %d):\n", hp_chain_energy(pt_best_chain(s)));
    hp_chain_print(pt_best_chain(s));
    printf("swap acceptance: %.4f\n", pt_swap_acceptance_rate(s));

    pt_free(s);
    return 0;
}
