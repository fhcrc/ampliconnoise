CC     = cc
CFLAGS = -O3
EFLAGS =  
EFILE  = FCluster
LIBS   = -lm 
OBJS   = FCluster.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : FCluster.c FCluster.h
	$(CC) $(CFLAGS) -c $*.c 

clean: 
	rm -rf *o FCluster
