CFLAGS        =  -g -I/usr/include/g++-2 -I../../mesh/shape -I../../analysis -I../../common -I../../external -I../../external/metis -I../../integral -I../../material -I../../matrix -I../../mesh -I../../multigrid -I../../post -I../../util -DPZENVIRONMENT

CC	      = g++

DEST	      = ${HOME}/bin

EXTHDRS	      = ../../common/pzreal.h \
		../../common/pzstream.h \
		../../material/pzmaterial.h \
		../../material/pzmathyperelastic.h \
		../../matrix/pzfmatrix.h \
		../../matrix/pzmatrix.h \
		../../util/pzvec.h

HDRS	      = checkconv.h \
		convhyper.h

LDFLAGS	      =  -L../../analysis -L../../integral -L../../material -L../../util -L../../external -L../../external/metis -L../../matrix -L../../mesh -L../../multigrid -L../../post -L../../matrix -L../../mesh/shape -lpzutil -lmultigrid -lanalysis -lexternal -lmesh -lshape -lpzutil -lquad -lmatrix -lmaterial -lpost -lmultigrid -lmetis -lmesh -lmatrix -lpthread 


LIBS	      =

LINKER	      = g++

MAKEFILE      = makefile

OBJS	      = convhyper.o \
		main.o

PRINT	      = pr

PROGRAM	      = geom

SRCS	      = convhyper.c \
		main.c

all:		$(PROGRAM)

dos2unix:;
		dos2unix *.c

$(PROGRAM):     $(OBJS) 
		$(LINKER) $(OBJS) $(LDFLAGS)  -o $(PROGRAM) 

clean:;		rm -f $(OBJS)

depend:;	mkmf -f $(MAKEFILE) PROGRAM=$(PROGRAM) DEST=$(DEST)

index:;		ctags -wx $(HDRS) $(SRCS)

install:	$(PROGRAM)
		install -s $(PROGRAM) $(DEST)

print:;		$(PRINT) $(HDRS) $(SRCS)

program:        $(PROGRAM)

tags:           $(HDRS) $(SRCS); ctags $(HDRS) $(SRCS)

update:		$(DEST)/$(PROGRAM)

$(DEST)/$(PROGRAM): $(SRCS) $(LIBS) $(HDRS) $(EXTHDRS)
		@make -f $(MAKEFILE) DEST=$(DEST) install
###
convhyper.o: convhyper.h ../../material/pzmathyperelastic.h \
	../../material/pzmaterial.h ../../common/pzreal.h \
	../../matrix/pzfmatrix.h ../../matrix/pzmatrix.h \
	../../common/pzstream.h ../../util/pzvec.h
main.o: convhyper.h ../../material/pzmathyperelastic.h \
	../../material/pzmaterial.h ../../common/pzreal.h \
	../../matrix/pzfmatrix.h ../../matrix/pzmatrix.h \
	../../common/pzstream.h ../../util/pzvec.h checkconv.h
