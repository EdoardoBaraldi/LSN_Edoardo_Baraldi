CC = g++
CFLAGS = -Wall -O3 --std=c++11

exe8.exe : exe8.o random.o
	$(CC) random.o exe8.o -o exe8.exe
exe8.o : exe8.cpp exe8.h
	$(CC) -c exe8.cpp -o exe8.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o exe8.exe
