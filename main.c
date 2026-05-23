#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "include/hp.h"
#include "include/mh.h"

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

static void print_usage(const char* prog) {
    fprintf(stderr,
        "Usage: %s [-s sequence] [-t temperature] [-n steps]\n"
        "\n"
        "  -s sequence   HP amino acid string (H/P), max length %d\n"
        "  -t temperature  Metropolis temperature (> 0), default 5.0\n"
        "  -n steps      Number of MH steps, default 50000\n"
        "  -h            Show this help\n"
        "\n"
        "Example:\n"
        "  %s -s HPPHPPHHPPHHPPHH -t 5.0 -n 50000\n",
        prog, HP_MAX_N, prog);
}

int main(int argc, char** argv) {
    const char* seq = "HPPHPPHHPPHHPPHH";
    double temperature = 5.0;
    int steps = 50000;

    int opt;
    while ((opt = getopt(argc, argv, "s:t:n:h")) != -1) {
        switch (opt) {
            case 's':
                seq = optarg;
                break;
            case 't':
                temperature = atof(optarg);
                break;
            case 'n':
                steps = atoi(optarg);
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
    if (temperature <= 0.0) {
        fprintf(stderr, "error: temperature must be > 0\n");
        return 1;
    }
    if (steps <= 0) {
        fprintf(stderr, "error: steps must be > 0\n");
        return 1;
    }

    hp_chain* chain = hp_chain_create(seq);
    if (chain == NULL) {
        fprintf(stderr, "error: hp_chain_create failed\n");
        return 1;
    }

    mh_sampler* s = mh_create(chain, temperature);
    if (s == NULL) {
        hp_chain_free(chain);
        fprintf(stderr, "error: mh_create failed\n");
        return 1;
    }

    printf("sequence:    %s\n", seq);
    printf("temperature: %.3f\n", temperature);
    printf("steps:       %d\n\n", steps);

    printf("Initial chain:\n");
    hp_chain_print(mh_chain(s));
    mh_run(s, steps);
    printf("Final chain:\n");
    hp_chain_print(mh_chain(s));
    printf("Best chain:\n");
    hp_chain_print(mh_best_chain(s));
    printf("best energy: %d\n", hp_chain_energy(mh_best_chain(s)));
    printf("acceptance:  %.2f\n", mh_acceptance_rate(s));

    mh_free(s);
    return 0;
}
