## makefile for tiling 

ODIR = ../obj
TDIR = ../builds
FILE = tiles
SFILE = $(FILE).c
SDIR = .
IRITDEF = -DLINUX386 -D__UNIX__ -DSTRICMP -DUSLEEP -DTIMES -DRAND -DIRIT_HAVE_XML_LIB

IRITLIB = -lIritExtLib  -lIritGrapLib -lIritUserLib -lIritRndrLib \
          -lIritBoolLib -lIritPrsrLib -lIritVMdlLib -lIritMdlLib \
          -lIritMvarLib -lIritTrimLib -lIritTrivLib -lIritTrngLib \
          -lIritSymbLib -lIritCagdLib -lIritGeomLib -lIritMiscLib \
          -lIritXtraLib
$(TDIR)/$(FILE):        $(ODIR)/$(FILE).o 
	$(CC)  -g -o $(TDIR)/$(FILE) $(ODIR)/$(FILE).o \
        	-L$(LDIR) -legads  -L$(IRIT_LIB) $(IRITLIB) $(IRITLIB) $(RPATH) -Wall -lm

$(ODIR)/$(FILE).o:	$(SDIR)/$(SFILE)  
	$(CC) -g -c $(COPTS) $(DEFINE) $(IRITDEF) -I$(IDIR)  -I. -I$(IRIT_INC) \
	$(SDIR)/$(SFILE) -o $(ODIR)/$(FILE).o

clean:





clean:
	(cd $(ODIR); rm -f $(FILE).o )

cleanall:	clean
	(cd $(LDIR); rm -f $(FILE).o )
