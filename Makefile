all: simulador

%:%.c
	cc -g -Wall -std=c99 -lm $< -o $@
