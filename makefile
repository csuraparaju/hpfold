CC      = gcc
CFLAGS  = -Wall -Wextra -Wno-unused-function -std=c99 -D_DEFAULT_SOURCE -I. -Iinclude
LDFLAGS = -lm

NVCC       = nvcc
CUDA_ARCH ?= native
NVCCFLAGS = -std=c++14 -O2 -Wno-deprecated-gpu-targets -arch=$(CUDA_ARCH) -I. -Iinclude

BUILD_DIR  = build

SRCS = src/hp.c src/moves.c src/mh.c src/rng.c main.c

C_SRCS_PT   = src/hp.c src/rng.c
CU_SRCS     = src/pt.cu
C_OBJS_PT   = $(addprefix $(BUILD_DIR)/, $(notdir $(C_SRCS_PT:.c=.o)))
CU_OBJS     = $(addprefix $(BUILD_DIR)/, $(notdir $(CU_SRCS:.cu=.o)))
MAIN_PT_OBJ = $(BUILD_DIR)/main_pt.o
PT_OBJS     = $(C_OBJS_PT) $(CU_OBJS) $(MAIN_PT_OBJ)

.PHONY: all mc_hp mc_hp_pt debug clean

all: mc_hp

mc_hp:
	$(CC) $(CFLAGS) -O2 $(SRCS) $(LDFLAGS) -o $@

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: src/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: src/%.cu | $(BUILD_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BUILD_DIR)/main_pt.o: main_pt.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

mc_hp_pt: $(PT_OBJS)
	$(NVCC) $(NVCCFLAGS) $^ $(LDFLAGS) -o $@

debug:
	$(CC) $(CFLAGS) -g -DDEBUG $(SRCS) $(LDFLAGS) -o mc_hp_dbg

clean:
	rm -f mc_hp mc_hp_dbg mc_hp_pt
	rm -rf $(BUILD_DIR)
	rm -f src/*.o main_pt.o
	rm -rf mc_hp.dSYM mc_hp_dbg.dSYM
