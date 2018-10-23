FC	= f77
FFLAGS	= -O 

INC	=
LDOPTS	=
CLIB	=
LIBS	=

TARGET	= libnumrec2.a
OBJS	= dqgaus.o dsplin2.o dsplie2.o dspline.o dsplint.o
      
$(TARGET): $(OBJS)
	ar r $@ $(OBJS)

.f.o:
	$(FC) $(FFLAGS) $(INC) -c $<

clean:
	@rm -f *.o *.a

.SUFFIXES: .o .f
