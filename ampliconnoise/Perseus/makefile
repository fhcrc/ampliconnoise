CC     = cc
CFLAGS = -O3 
EFLAGS =  
EFILE  = Perseus
LIBS   = -lm -lgsl -lgslcblas
OBJS   = Perseus.o 

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : Perseus.c Perseus.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o Perseus
