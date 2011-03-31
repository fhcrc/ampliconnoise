CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = SeqNoise
LIBS   = -lm 
OBJS   = SeqNoise.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SeqNoise.c SeqNoise.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SeqNoise
