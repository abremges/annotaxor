CC=gcc
CFLAGS=-Wall -pedantic -std=gnu99 -O3

annotaxor:annotaxor.c kseq.h khash.h
	$(CC) $(CFLAGS) annotaxor.c -o $@ -lz

clean:
	rm -fr gmon.out *.o ext/*.o a.out annotaxor *~ *.a *.dSYM session*
