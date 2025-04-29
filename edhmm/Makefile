CC = gcc

# debug for memory CFLAGS = -fsanitize=address -g
CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o

all: edhmm

edhmm: $(OBJS)
	$(CC) $(CFLAGS) -o edhmm $(OBJS) -lm
	@chmod +x edhmm

main.o: main.c parser.c model.c model.h
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compiling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compiling model"
	$(CC) $(CFLAGS) -c model.c

clean:
	@echo "Removing everything but the source file"
	rm -f edhmm *.o

# with default debug 
debug: edhmm
	lldb ./edhmm seq

# with default seq input
test: edhmm
	./edhmm seq
	
.PHONY: all clean debug