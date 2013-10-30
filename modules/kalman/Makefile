CC     = clang
CFLAGS = -Wall
CFILES = main.c kalman_matrix.c
OFILES = $(CFILES:.c=.o)
EXE = kalman_matrix_test.x

LIBS   =  


$(EXE): $(OFILES)
	$(CC) $(OFILES) ${LIBS} -o $(EXE)


main.o: main.c
	$(CC) -c main.c $(CFLAGS)

kalman_matrix.o: kalman_matrix.c
	$(CC) -c kalman_matrix.c $(CFLAGS)


debug:
	$(CC) $(CFILES) ${LIBS} -g -o debug.x


clean:
	@/bin/rm -f *.o *.x *.c~ *.h~


main.o: main.c kalman_matrix.h
kalman_matrix.o: kalman_matrix.c kalman_matrix.h
