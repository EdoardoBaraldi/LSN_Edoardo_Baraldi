CC = /home/edoardo/MPICH/bin/mpicxx
CFLAGS = -Wall -O3 --std=c++17

exe10.exe : exe10.o random.o genetic_alg.o
	$(CC) random.o genetic_alg.o exe10.o -o exe10.exe
exe10.o : exe10.cpp exe10.h
	$(CC) -c exe10.cpp -o exe10.o $(CFLAGS)
genetic_alg.o : genetic_alg.cpp genetic_alg.h 
	$(CC) -c genetic_alg.cpp -o genetic_alg.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o exe10.exe
