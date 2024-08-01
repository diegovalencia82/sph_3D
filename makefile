GFORTFLAGS = -fno-second-underscore

LIBS = -L$(NEMOLIB)

clean:
	rm -f *.o *.a core 

.f:

#	gfortran $(GFORTFLAGS) -o $* $*.f -L/usr/local/pgplot -L/usr/X11/lib -lpgplot -lX11 ./mysph.a

	gfortran $(GFORTFLAGS) -o $* $*.f ./mysph.a

	rm -f *.o









#gfortran -o join_two_models join_two_models.f -L/usr/local/pgplot -L/usr/X11/lib -lpgplot -lX11
