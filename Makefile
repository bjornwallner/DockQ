FNAT=molecule.o fnat.o

#CCFLAG=-static -static-libgcc 
#CCFLAG=-w 
LFLAG=-O3 -funroll-loops -Isrc/


fnat: $(FNAT)
	$(CC) $(LFLAG) $(CCFLAG) -o fnat $(FNAT) -lm 

fnat.o: src/fnat.c
	$(CC) $(LFLAG) $(CCFLAG) -c src/fnat.c -lm 
molecule.o: src/molecule.c src/molecule.h 
	$(CC) $(LFLAG) $(CCFLAG) -c src/molecule.c -lm 

.c.o:	
	$(CC) -c $(LFLAG) $(CCFLAG)  $*.c -lm 


