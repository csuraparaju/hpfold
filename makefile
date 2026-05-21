make:
	gcc -Wall -Wextra  -std=c99 -O2 hp.c main.c -lm -o mc_hp
clean:
	rm mc_hp