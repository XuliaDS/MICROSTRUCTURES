#
IDIR = $(ESP_ROOT)/include
include $(IDIR)/$(ESP_ARCH)
LDIR = $(ESP_ROOT)/lib
#ifdef ESP_BLOC
ODIR = .
TDIR = .
#else
#ODIR = .
#TDIR = $(ESP_ROOT)/bin
#endif
FILE  = trimmingApproximation
SFILE = $(FILE).c
SDIR = .

IRITDEF = -DLINUX386 -D__UNIX__ -DSTRICMP -DUSLEEP -DTIMES -DRAND -DIRIT_HAVE_XML_LIB

IRITLIB = -lIritExtLib  -lIritGrapLib -lIritUserLib -lIritRndrLib \
          -lIritBoolLib -lIritPrsrLib -lIritVMdlLib -lIritMdlLib \
          -lIritMvarLib -lIritTrimLib -lIritTrivLib -lIritTrngLib \
          -lIritSymbLib -lIritCagdLib -lIritGeomLib -lIritMiscLib \
          -lIritXtraLib

$(TDIR)/$(FILE):	$(ODIR)/$(FILE).o  
	$(CC) -g -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
		-L$(LDIR)  -legads   -L$(IRIT_LIB) $(IRITLIB) $(IRITLIB) -Wall -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE) $(IDIR)/egads.h $(IDIR)/egadsTypes.h \
			$(IDIR)/egadsErrors.h 
	$(CC)  -g -c $(COPTS) $(DEFINE)  $(IRITDEF) -I$(IDIR)  -I. -I$(IRIT_INC) \
	$(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:
	-rm $(ODIR)/$(FILE).o

cleanall:	clean
	-rm $(TDIR)/$(FILE)
