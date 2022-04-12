F77 = gfortran-mp-11
FFLAGS = -c -O2 -std=legacy
LFLAGS = -O2 -Wall -Wsurprising -fbounds-check -std=legacy
HOSTLIBS = -L/opt/local/lib/gcc11 -lm -lgcc  
FITLIB = -L./cfitsio/lib -lcfitsio -lcurl 
FSOFLAGS = -shared -O2 -std=legacy -Wall -Wsurprising -fbounds-check


# RULES:
.SUFFIXES: .o .f
.f.o:
	$(F77) $(FFLAGS) $<

daophot: daophot.o pckpsf.o find.o fotometry.o \
        psf.o peak.o nstar-mcmc.o fudge.o addstar.o substar.o \
        group.o sort.o lnxsubs.o fitsubs.o iosubs.o mathsubs.o
	$(F77) $(LFLAGS) -o daophot daophot.o pckpsf.o find.o fotometry.o \
        psf.o peak.o nstar-mcmc.o fudge.o addstar.o substar.o \
        group.o sort.o lnxsubs.o fitsubs.o \
	iosubs.o mathsubs.o \
        $(HOSTLIBS) $(FITLIB)

allstar: allstar.o allstsubs.o lnxsubs.o \
	iosubs.o mathsubs.o fitsubs.o
	$(F77) $(LFLAGS) -o allstar allstar.o allstsubs.o \
	lnxsubs.o iosubs.o mathsubs.o fitsubs.o \
	$(HOSTLIBS) $(FITLIB)

daomaster: daomaster.o iosubs.o mathsubs.o lnxsubs.o dummysm.o
	$(F77) $(LFLAGS) -o daomaster daomaster.o iosubs.o mathsubs.o lnxsubs.o dummysm.o

daomatch: daomatch.o iosubs.o mathsubs.o lnxsubs.o
	$(F77) $(LFLAGS) -o daomatch daomatch.o iosubs.o mathsubs.o lnxsubs.o

montage2: montage2.o mathsubs.o iosubs.o lnxsubs.o fitsubs.o
	$(F77) $(LFLAGS) -o montage2 montage2.o mathsubs.o iosubs.o \
	lnxsubs.o fitsubs.o \
        $(HOSTLIBS) $(FITLIB)

allframe: allframe.o fitsubs.o lnxsubs.o iosubs.o mathsubs.o
	$(F77) $(LFLAGS) -o allframe allframe.o fitsubs.o \
	lnxsubs.o iosubs.o mathsubs.o \
        $(HOSTLIBS) $(FITLIB)

daogrow: daogrow.o iosubs.o mathsubs.o lnxsubs.o
	$(F77) $(LFLAGS) -o daogrow daogrow.o iosubs.o mathsubs.o lnxsubs.o

pywrapper: daophot.o pckpsf.o find.o fotometry.o \
        psf.o peak.o nstar-mcmc_pywrapper.o fudge.o addstar.o substar.o \
        group.o sort.o lnxsubs.o fitsubs.o iosubs.o mathsubs.o
	$(F77) $(FSOFLAGS) $(HOSTLIBS) $(FITLIB) pywrapper.f -o pywrapper.so \
	daophot.o pckpsf.o find.o fotometry.o \
        psf.o peak.o nstar-mcmc_pywrapper.o fudge.o addstar.o substar.o \
        group.o sort.o lnxsubs.o fitsubs.o \
	iosubs.o mathsubs.o
