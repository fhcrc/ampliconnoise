CC     = mpicc
CFLAGS = -O3
EFLAGS = 
EFILE  = PyroNoiseM
LIBS   =  -lm
OBJS   = PyroNoiseM.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : PyroNoiseM.c PyroNoiseM.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o PyroNoiseM
