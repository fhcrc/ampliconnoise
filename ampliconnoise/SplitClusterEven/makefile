CC     = cc
CFLAGS = -O3
EFLAGS =  
EFILE  = SplitClusterEven
LIBS   = -lm 
OBJS   = SplitClusterEven.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SplitClusterEven.c SplitClusterEven.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SplitClusterEven
