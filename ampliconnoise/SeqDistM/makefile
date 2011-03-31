CC     = mpicc
CFLAGS = -O3
EFLAGS =  
EFILE  = SeqDistM
LIBS   = -lm 
OBJS   = SeqDistM.o

$(EFILE) : $(OBJS)
	@echo "linking..."
	$(CC) $(EFLAGS) -o $(EFILE) $(OBJS) $(LIBS)

$(OBJS) : SeqDistM.c SeqDistM.h
	$(CC) $(CFLAGS) -c $*.c 

clean:
	rm -rf *.o SeqDistM
