##

CC = g++
OPT = -std=c++11 -static-libstdc++
LDFLAGS = -lm

HDRS = utils.hpp md_vector.h
SRCS = main.cpp utils.cpp
OBJS = main.o utils.o

PROGRAM = naiv_seq

##  link
$(PROGRAM) : $(OBJS)
	$(CC) -o $(PROGRAM) $(OBJS) $(OPT) $(LDFLAGS)

##  dependency
main.o : main.cpp utils.hpp md_vector.h
	$(CC) -c main.cpp $(OPT)

utils.o : utils.cpp utils.hpp
	$(CC) -c utils.cpp $(OPT)

.PHONY: clean
clean:
	rm -f *.o a.out core $(PROGRAM)

.PHONY: all
all: clean $(PROGRAM)
