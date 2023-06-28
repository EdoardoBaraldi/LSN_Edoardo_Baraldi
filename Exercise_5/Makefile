CC = g++
CFLAGS = -Wall -O3 --std=c++11

exe5.exe : exe5.o random.o
	$(CC) random.o exe5.o -o exe5.exe  	
exe5.o : exe5.cpp exe5.h
	$(CC) -c exe5.cpp -o exe5.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
