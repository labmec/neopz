CFLAGS        =  -g -I/usr/include/g++-2 -I../analysis -I../common -I../external -I../external/metis -I../integral -I../material -I../matrix -I../mesh -I../multigrid -I../post -I../util -DPZENVIRONMENT

CC	      = g++

DEST	      = ${HOME}/bin

EXTHDRS	      = ../analysis/pzanalysis.h \
		../common/pzerror.h \
		../common/pzreal.h \
		../common/pzstream.h \
		../external/pzmetis.h \
		../external/pzrenumbering.h \
		../integral/pzquad.h \
		../material/pzbndcond.h \
		../material/pzelasmat.h \
		../material/pzmat1dlin.h \
		../material/pzmat2dlin.h \
		../material/pzmaterial.h \
		../material/pzmatplaca2.h \
		../material/pzmultplaca.h \
		../material/tpzmultcamada.h \
		../matrix/pzblock.h \
		../matrix/pzfmatrix.h \
		../matrix/pzmatrix.h \
		../matrix/pzskylmat.h \
		../matrix/pzsolve.h \
		../mesh/pzcmesh.h \
		../mesh/pzcompel.h \
		../mesh/pzconnect.h \
		../mesh/pzelct2d.h \
		../mesh/pzelg1d.h \
		../mesh/pzelgq2d.h \
		../mesh/pzelgt2d.h \
		../mesh/pzelmat.h \
		../mesh/pzgeoel.h \
		../mesh/pzgmesh.h \
		../mesh/pzgnode.h \
		../mesh/pzintel.h \
		../multigrid/pztrnsform.h \
		../util/pzadmchunk.h \
		../util/pzchunk.h \
		../util/pzmanvector.h \
		../util/pzstack.h \
		../util/pzvec.h

HDRS	      =

LDFLAGS	      = -L../analysis -L../integral -L../material -L../util -L../external -L../external/metis -L../matrix -L../mesh -L../multigrid -L../post -L../matrix -L../mesh/shape -lmultigrid -lpost -lanalysis -lexternal -lmesh -lshape -lpzutil -lquad -lmatrix -lmaterial -lpost -lmultigrid -lmetis -lmatrix


LIBS	      = ../multigrid/libmultigrid.a \
		../analysis/libanalysis.a \
		../external/metis/libmetis.a \
		../external/libexternal.a \
		../mesh/libmesh.a \
		../util/libpzutil.a \
		../integral/libquad.a \
		../post/libpost.a \
		../matrix/libmatrix.a \
		../material/libmaterial.a

LINKER	      = g++

MAKEFILE      = makefile

OBJS	      = jornadas.o

PRINT	      = pr

PROGRAM	      = jornadas

SRCS	      = jornadas.c

all:		$(PROGRAM)

dos2unix:;
		dos2unix *.c

$(PROGRAM):     $(OBJS) 
		$(LINKER) $(OBJS) $(LDFLAGS)  -o $(PROGRAM) -lpzutil

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
jornadas.o: ../matrix/pzfmatrix.h ../matrix/pzmatrix.h ../common/pzstream.h \
	../common/pzreal.h ../matrix/pzskylmat.h ../util/pzvec.h \
	../util/pzmanvector.h ../mesh/pzgmesh.h ../mesh/pzgnode.h \
	../material/pzbndcond.h ../material/pzmaterial.h ../mesh/pzcmesh.h \
	../util/pzadmchunk.h ../util/pzchunk.h ../util/pzstack.h \
	../matrix/pzblock.h ../mesh/pzconnect.h ../analysis/pzanalysis.h \
	../matrix/pzsolve.h ../mesh/pzelg1d.h ../mesh/pzgeoel.h \
	../common/pzerror.h ../multigrid/pztrnsform.h \
	../material/pzmat1dlin.h ../mesh/pzelmat.h ../material/pzelasmat.h \
	../mesh/pzcompel.h ../mesh/pzelgq2d.h ../mesh/pzelgt2d.h \
	../material/pzmat2dlin.h ../external/pzmetis.h \
	../external/pzrenumbering.h ../material/pzmatplaca2.h \
	../material/pzmultplaca.h ../material/tpzmultcamada.h \
	../mesh/pzelct2d.h ../mesh/pzintel.h ../integral/pzquad.h
