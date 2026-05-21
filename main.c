#include <stdio.h>
#include <string.h>
#include "hp.h"

static void run_test(const char* seq, int steps, double temp) {
    printf("Sequence: %s\n", seq);
    printf("Length:   %zu\n", strlen(seq));
    printf("Temp:     %.3f\n", temp);
    printf("Steps:    %d\n", steps);

    mc_trajectory* traj = mc_traj_create(seq);

    printf("\nInitial configuration:\n");
    mc_traj_print(traj);

    printf("Initial energy: %d\n", mc_traj_energy(traj));

    int best = mc_traj_run(traj, steps, temp);

    printf("\nFinal configuration:\n");
    mc_traj_print(traj);

    printf("Final energy: %d\n", mc_traj_energy(traj));

    printf("Best energy encountered: %d\n", best);
    printf("\n");
    printf("\n");

    mc_traj_free(traj);

}


int main(void) {
    run_test(
        "HPPHPPHHPPHHPPHH",
        50000,
        5.0
    );

    return 0;

}