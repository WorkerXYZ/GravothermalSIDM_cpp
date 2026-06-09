CC ?= gcc
EXE := gravothermal_sidm
SRC := gravothermal_sidm.c
HDR := gravothermal_sidm.h

USE_OPENMP ?= 1
OMP_MIN_SHELLS ?= 2048
ARGS ?=

CFLAGS ?= -O3 -std=c11 -Wall -Wextra -pedantic
LDFLAGS ?=

ifeq ($(USE_OPENMP),1)
	CFLAGS += -fopenmp
	LDFLAGS += -fopenmp
endif

CFLAGS += -DGSIDM_OPENMP_MIN_SHELLS=$(OMP_MIN_SHELLS)
LDFLAGS += -lm

.DEFAULT_GOAL := build

.PHONY: all build run clean

all: build

build: $(EXE)

$(EXE): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $(SRC) $(LDFLAGS)

run: $(EXE)
	./$(EXE) $(ARGS)

clean:
	rm -f $(EXE)
