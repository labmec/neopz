CFLAGS        =  -g -I/usr/include/g++-2 -I../../mesh/shape -I../../analysis -I../../common -I../../external -I../../external/metis -I../../integral -I../../material -I../../matrix -I../../mesh -I../../multigrid -I../../post -I../../util -DPZENVIRONMENT

CC	      = g++

DEST	      = ${HOME}/bin

EXTHDRS	      = ../../analysis/pzanalysis.h \
		../../common/pzerror.h \
		../../common/pzreal.h \
		../../common/pzstream.h \
		../../material/pzbndcond.h \
		../../material/pzmat2dlin.h \
		../../material/pzmaterial.h \
		../../matrix/pzblock.h \
		../../matrix/pzfmatrix.h \
		../../matrix/pzmatrix.h \
		../../mesh/pzcmesh.h \
		../../mesh/pzcompel.h \
		../../mesh/pzconnect.h \
		../../mesh/pzelgq2d.h \
		../../mesh/pzgeoel.h \
		../../mesh/pzgmesh.h \
		../../mesh/pzgnode.h \
		../../mesh/shape/pzshapelinear.h \
		../../multigrid/pztrnsform.h \
		../../post/pzdxmesh.h \
		../../post/pzgraphmesh.h \
		../../util/pzadmchunk.h \
		../../util/pzchunk.h \
		../../util/pzmanvector.h \
		../../util/pzstack.h \
		../../util/pzvec.h

HDRS	      =

LDFLAGS	      =  -L../../analysis -L../../integral -L../../material -L../../util -L../../external -L../../external/metis -L../../matrix -L../../mesh -L../../multigrid -L../../post -L../../matrix -L../../mesh/shape -lpzutil -lmultigrid -lanalysis -lexternal -lmesh -lshape -lpzutil -lquad -lmatrix -lmaterial -lpost -lmultigrid -lmetis -lmesh -lmatrix -lpthread 


LIBS	      =

LINKER	      = g++

MAKEFILE      = Makefile

OBJS	      = geom.o

PRINT	      = pr

PROGRAM	      = geom

SRCS	      = geom.c

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
geom.o: ../../mesh/pzgmesh.h ../../common/pzreal.h ../../mesh/pzgnode.h \
	../../material/pzbndcond.h ../../material/pzmaterial.h \
	../../util/pzvec.h ../../util/pzmanvector.h ../../matrix/pzfmatrix.h \
	../../matrix/pzmatrix.h ../../common/pzstream.h ../../mesh/pzcmesh.h \
	../../util/pzadmchunk.h ../../util/pzchunk.h ../../util/pzstack.h \
	../../matrix/pzblock.h ../../mesh/pzconnect.h \
	../../analysis/pzanalysis.h ../../mesh/pzelgq2d.h \
	../../mesh/pzgeoel.h ../../common/pzerror.h \
	../../multigrid/pztrnsform.h ../../material/pzmat2dlin.h \
	../../post/pzdxmesh.h ../../post/pzgraphmesh.h ../../mesh/pzcompel.h \
	../../mesh/shape/pzshapelinear.h
