CFLAGS        =  -g -I/usr/include/g++-2 -I../../mesh/shape -I../../analysis -I../../common -I../../external -I../../external/metis -I../../integral -I../../material -I../../matrix -I../../mesh -I../../multigrid -I../../post -I../../util -DPZENVIRONMENT

CC	      = g++

DEST	      = ${HOME}/bin

EXTHDRS	      = ../../analysis/pzanalysis.h \
		../../analysis/pznonlinanalysis.h \
		../../common/pzerror.h \
		../../common/pzreal.h \
		../../common/pzstream.h \
		../../integral/pzquad.h \
		../../material/pzbndcond.h \
		../../material/pzelasmat.h \
		../../material/pzmat1dlin.h \
		../../material/pzmat2dlin.h \
		../../material/pzmaterial.h \
		../../material/pzmathyperelastic.h \
		../../material/pzmattest.h \
		../../material/pzmattest3d.h \
		../../matrix/pzblock.h \
		../../matrix/pzfmatrix.h \
		../../matrix/pzmatrix.h \
		../../matrix/pzskylmat.h \
		../../matrix/pzsolve.h \
		../../matrix/pztempmat.h \
		../../mesh/pzcmesh.h \
		../../mesh/pzcompel.h \
		../../mesh/pzconnect.h \
		../../mesh/pzelc1d.h \
		../../mesh/pzelcc3d.h \
		../../mesh/pzelcpi3d.h \
		../../mesh/pzelcpr3d.h \
		../../mesh/pzelct2d.h \
		../../mesh/pzelct3d.h \
		../../mesh/pzelg1d.h \
		../../mesh/pzelgc3d.h \
		../../mesh/pzelgpi3d.h \
		../../mesh/pzelgpr3d.h \
		../../mesh/pzelgq2d.h \
		../../mesh/pzelgt2d.h \
		../../mesh/pzelgt3d.h \
		../../mesh/pzelmat.h \
		../../mesh/pzgeoel.h \
		../../mesh/pzgmesh.h \
		../../mesh/pzgnode.h \
		../../mesh/pzintel.h \
		../../mesh/shape/pzshapelinear.h \
		../../multigrid/pztrnsform.h \
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

OBJS	      = hyperelastic.o

PRINT	      = pr

PROGRAM	      = testmat

SRCS	      = hyperelastic.c

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
hyperelastic.o: ../../mesh/pzgeoel.h ../../common/pzerror.h \
	../../common/pzreal.h ../../mesh/pzgmesh.h ../../mesh/pzgnode.h \
	../../material/pzbndcond.h ../../material/pzmaterial.h \
	../../util/pzvec.h ../../util/pzmanvector.h ../../matrix/pzfmatrix.h \
	../../matrix/pzmatrix.h ../../common/pzstream.h ../../mesh/pzcmesh.h \
	../../util/pzadmchunk.h ../../util/pzchunk.h ../../util/pzstack.h \
	../../matrix/pzblock.h ../../mesh/pzconnect.h \
	../../analysis/pzanalysis.h ../../multigrid/pztrnsform.h \
	../../matrix/pzskylmat.h ../../analysis/pznonlinanalysis.h \
	../../mesh/pzcompel.h ../../mesh/shape/pzshapelinear.h \
	../../matrix/pztempmat.h ../../matrix/pzsolve.h ../../mesh/pzintel.h \
	../../mesh/pzelg1d.h ../../mesh/pzelc1d.h ../../integral/pzquad.h \
	../../mesh/pzelgq2d.h ../../mesh/pzelgt2d.h ../../mesh/pzelct2d.h \
	../../mesh/pzelgc3d.h ../../mesh/pzelcc3d.h ../../mesh/pzelgt3d.h \
	../../mesh/pzelct3d.h ../../mesh/pzelgpi3d.h ../../mesh/pzelcpi3d.h \
	../../mesh/pzelgpr3d.h ../../mesh/pzelcpr3d.h ../../mesh/pzelmat.h \
	../../material/pzelasmat.h ../../material/pzmattest.h \
	../../material/pzmat1dlin.h ../../material/pzmat2dlin.h \
	../../material/pzmattest3d.h ../../material/pzmathyperelastic.h
