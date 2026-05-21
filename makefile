all:
	gcc -Wall -Wextra -std=c99 -O2 hp.c main.c -lm -o mc_hp

debug:
	gcc -Wall -Wextra -DDEBUG -std=c99 -g hp.c main.c -lm -o mc_hp_dbg

clean:
	rm -rf mc_hp_dbg.dSYM
	rm -rf mc_hp.dSYM
	rm -f mc_hp mc_hp_dbg