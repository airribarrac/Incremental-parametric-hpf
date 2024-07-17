CFLAGS=-O4 -g  -DFAST_IO -DDIRECTED_CASE 
BINDIR=bin
SOURCES = src/1.0/incremental.c
OBJECTS = $(SOURCES:.c=.o)
TARGET = incremental
CC=gcc

all: $(TARGET) 

clean:
	rm -f $(OBJECTS) $(TARGET) 

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o bin/$(TARGET) $(OBJECTS)


$(TARGET)_anytime: $(OBJECTS)
	$(CC) $(CFLAGS) -DANYTIME -o bin/$(TARGET)_anytime $(SOURCES)
