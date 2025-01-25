#!/usr/bin/env bash

# PoC Only

nasm -O2 -o dft_sse2.o -f elf64 dft.asm
gcc -O3 -o test main.c dft_sse2.o -lm
