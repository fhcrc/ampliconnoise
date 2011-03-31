CC     = cc
CFLAGS = -O3
EFLAGS =  
EFILE  = FClusterM
LIBS   = -lm 
OBJS   = FClusterM.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : FClusterM.c FClusterM.h
	$(CC) $(CFLAGS) -c $*.c 

clean: 
	rm -rf *o FClusterM
