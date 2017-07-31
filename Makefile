CFLAGS=-g -Wall
#CFLAGS=-O2 -Wall


run: graphecycle
	./graphecycle  4672 0 20 600 1

val: graphecycle
	valgrind --track-origins=yes ./graphecycle  4672 0 20 600 1

graphecycle: graphecycle.o lecture_molecule_sdf.o fonctions_molecules.o 
	gcc ${CFLAGS} graphecycle.o lecture_molecule_sdf.o fonctions_molecules.o -o graphecycle

graphecycle.o: graphecycle.c graphecycle.h
	gcc ${CFLAGS} -c graphecycle.c

fonctions_molecules.o: fonctions_molecules.c graphecycle.h
		gcc ${CFLAGS} -c fonctions_molecules.c
		
	
lecture_molecule_sdf.o: lecture_molecule_sdf.c lecture_molecule_sdf.h
	gcc ${CFLAGS} -c lecture_molecule_sdf.c

clean: 
	rm -f graphecycle
	rm -f graphecycle.o
	rm -f lecture_molecule_sdf.o
	rm -f fonctions_molecules.o

