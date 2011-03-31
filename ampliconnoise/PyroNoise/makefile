CC     = mpicc
CFLAGS = -O3
EFLAGS = 
EFILE  = PyroNoise
LIBS   =  -lm
OBJS   = PyroNoise.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : PyroNoise.c PyroNoise.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o PyroNoise
