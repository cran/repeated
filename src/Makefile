LIBNAME=repeated
LD=ld

OBJS= chidden.o cutil.o dist.o gar.o hidden.o inthaz.o kcountb.o \
	kserieb.o logitord.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(LD) $(SHLIBLDFLAGS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
