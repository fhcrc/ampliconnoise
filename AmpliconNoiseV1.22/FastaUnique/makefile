CC     = cc
CFLAGS = -O3 
EFLAGS =  
EFILE  = FastaUnique
LIBS   = -lm 
OBJS   = FastaUnique.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : FastaUnique.c FastaUnique.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o FastaUnique
