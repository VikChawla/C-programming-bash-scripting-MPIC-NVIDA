CC=gcc
<<<<<<< HEAD
CFLAGS=-lm -O3
=======
CFLAGS=-lm
>>>>>>> 3dee7bb1001e17039d938df952127c7b4574a7a3
DEPS = intersections.h
OBJ = intersections.o \
      utility.o \
      main.o 

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

intersections: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f ./*.o ./*~ intersections
