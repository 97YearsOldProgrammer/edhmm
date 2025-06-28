CC = gcc

#CFLAGS = -fsanitize=address -g
CFLAGS = -g -Wall -Wextra	# normal with default debug; mostly useless for runtime bug

OBJS = main.o parser.o model.o output.o

all: edhmm2

edhmm2: $(OBJS)
	$(CC) $(CFLAGS) -o edhmm2 $(OBJS) -lm
	@chmod +x edhmm2

main.o: main.c parser.c model.c model.h output.c
	@echo "Compiling main"
	$(CC) $(CFLAGS) -c main.c

parser.o: parser.c model.h
	@echo "Compiling function and data structure for computation"
	$(CC) $(CFLAGS) -c parser.c

model.o: model.c model.h
	@echo "Compiling model"
	$(CC) $(CFLAGS) -c model.c

output.o: output.c model.h
	@echo "Compiling output session"
	$(CC) $(CFLAGS) -c output.c

clean:
	@echo "Removing everything but the source file"
	rm -f edhmm2 *.o

# with default debug 
debug: edhmm2
	lldb ./edhmm2 seq

# with default seq input
test: edhmm2
	./edhmm2 seq 1 ../models/acc.pwm ../models/don.pwm ../models/exon.mm ../models/intron.mm ../models/exon.len ../models/intron.len
	
.PHONY: all clean debug