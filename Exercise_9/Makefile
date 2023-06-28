CC = g++
CFLAGS = -Wall -O3 --std=c++17

exe9.exe : exe9.o random.o genetic_alg.o
	$(CC) random.o genetic_alg.o exe9.o -o exe9.exe
exe9.o : exe9.cpp exe9.h
	$(CC) -c exe9.cpp -o exe9.o $(CFLAGS)
genetic_alg.o : genetic_alg.cpp genetic_alg.h 
	$(CC) -c genetic_alg.cpp -o genetic_alg.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o exe9.exe
