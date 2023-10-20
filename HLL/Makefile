CC = gcc
CFLAGS = -Wall -Werror
LIBS = -lm

SRCS = intMain.c EOS_cold.c read_par_ic.c read_par_eos.c auxFun.c highOrder.c itmethod.c HLLrel.c reconstructions.c rel_euler.c
OBJS = $(SRCS:.c=.o)
EXECUTABLE = Rsolver_HLL

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJS)
