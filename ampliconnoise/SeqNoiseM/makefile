CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = SeqNoiseM
LIBS   = -lm 
OBJS   = SeqNoiseM.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SeqNoiseM.c SeqNoiseM.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SeqNoise
