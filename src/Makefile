LIBNAME=repeated

OBJS= chidden.o cutil.o dist.o gar.o hidden.o inthaz.o kcountb.o \
	kserieb.o logitord.o binnest.o

$(LIBNAME).$(SHLIBEXT):	$(OBJS)
	$(SHLIBLD) $(SHLIBLDFLAGS) $(FLIBS) -o $@ $(OBJS)

clean:
	@rm -f *.o *.$(SHLIBEXT)

realclean: clean
