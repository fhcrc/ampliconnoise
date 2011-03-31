CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = SeqDist
LIBS   = -lm 
OBJS   = SeqDist.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SeqDist.c SeqDist.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SeqDist
