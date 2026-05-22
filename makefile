CC      = gcc
CFLAGS  = -Wall -Wextra -std=c99 -I. -Iinclude
LDFLAGS = -lm

SRCS = src/hp.c src/moves.c src/mh.c src/rng.c main.c

.PHONY: all debug clean

all: mc_hp

mc_hp:
	$(CC) $(CFLAGS) -O2 $(SRCS) $(LDFLAGS) -o $@

debug:
	$(CC) $(CFLAGS) -g -DDEBUG $(SRCS) $(LDFLAGS) -o mc_hp_dbg

clean:
	rm -f mc_hp mc_hp_dbg
	rm -rf mc_hp.dSYM mc_hp_dbg.dSYM