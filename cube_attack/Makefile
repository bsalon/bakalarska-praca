CC=gcc
CFLAGS += -std=c18 -g -Wall -Wextra

all: degree_3 degree_5

bool_degree_3.o: bool_degree_3.c
	$(CC) $(CFLAGS) -c bool_degree_3.c

degree_3: bool_degree_3.o
	$(CC) -o degree_3 bool_degree_3.o

bool_degree_5.o: bool_degree_5.c
	$(CC) $(CFLAGS) -c bool_degree_5.c

degree_5: bool_degree_5.o
	$(CC) -o degree_5 bool_degree_5.o

.PHONY: clean

clean:
	rm -f degree_3 degree_5 bool_degree_3.o bool_degree_5.o
