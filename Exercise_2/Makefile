CC = g++
CFLAGS = -Wall -O3 --std=c++11

exe1.exe : exe1.o random.o
	$(CC) random.o exe1.o -o exe1.exe
exe2.exe : exe2.o random.o
	$(CC) random.o exe2.o -o exe2.exe  	
exe1.o : exe1.cpp mylib.h
	$(CC) -c exe1.cpp -o exe1.o $(CFLAGS)
exe2.o : exe2.cpp mylib.h
	$(CC) -c exe2.cpp -o exe2.o $(CFLAGS)	
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o *.exe *.out
