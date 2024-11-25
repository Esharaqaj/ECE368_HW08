CC = gcc
CFLAGS = -Wall -Wextra -O2
TARGET = a8

all: $(TARGET)

$(TARGET): main.o
	$(CC) $(CFLAGS) -o $(TARGET) main.o

main.o: main.c main.h
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f *.o $(TARGET)
