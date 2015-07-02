CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = PyroDist
LIBS   = -lm 
OBJS   = PyroDist.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : PyroDist.c PyroDist.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o PyroDist
