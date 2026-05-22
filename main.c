#include <stdio.h>
#include <string.h>
#include "include/mh.h"


int main(void) {
  mh_sampler* s = mh_create(hp_chain_create("HPPHPPHHPPHHPPHH"), 5.0);
  printf("Initial chain:\n");
  hp_chain_print(mh_chain(s));           // initial
  mh_run(s, 50000);
  printf("Final chain:\n");
  hp_chain_print(mh_chain(s));           // final
  printf("Best chain:\n");
  hp_chain_print(mh_best_chain(s));      // best seen
  printf("best energy: %d\n", hp_chain_energy(mh_best_chain(s)));
  printf("acceptance:  %.2f\n", mh_acceptance_rate(s));
  mh_free(s);

}