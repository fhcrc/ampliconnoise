CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = NDist
LIBS   = -lm  
OBJS   = NDist.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : NDist.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o NDist
