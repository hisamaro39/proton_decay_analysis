CC	  = gcc
CFLAGS	  =  -g -o
LOADLIBES = -lm
#LIBDIR    = /u/nurock/kenkou/work/sk/ndecay/gsl
#LIBDIR = /home/jsjang/gsl
#LIBDIR = /net/sukatmd1/work21/nishino/gsl
LIBDIR = /usr/include/gsl
INCGSL    = -I$(LIBDIR)/include/gsl -I$(LIBDIR)/include
LIBGSL    = -L$(LIBDIR)/lib -lgsl -lgslcblas
#OBJ       = poisson_cl.o 

calclimit: calclimit.o $(OBJ) 
	$(CC) $(CFLAGS) ./calclimit calclimit.o $(OBJ) $(LIBGSL) $(LOADLIBES)

make_table: make_table.o $(OBJ) 
	$(CC) $(CFLAGS) ./make_table make_table.o $(OBJ) $(LIBGSL) $(LOADLIBES)

	
calclimit_sk12: calclimit_sk12.o $(OBJ) 
	$(CC) $(CFLAGS) ./calclimit_sk12 calclimit_sk12.o $(OBJ) $(LIBGSL) $(LOADLIBES)

combinelimit: combinelimit.o $(OBJ) 
	$(CC) $(CFLAGS) ./combinelimit combinelimit.o $(OBJ) $(LIBGSL) $(LOADLIBES)

.c.o:
	$(CC) $(INCGSL) -c $<
clean:
	@rm -f *.o *~ core *.trace fort.* *.log

