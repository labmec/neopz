cmake_minimum_required (VERSION 2.6)
project (PZ)

include_directories(
	"${PROJECT_SOURCE_DIR}/Util"
	"${PROJECT_SOURCE_DIR}/Common"
	"${PROJECT_SOURCE_DIR}/Save"
	"${PROJECT_SOURCE_DIR}/Matrix"
	"${PROJECT_SOURCE_DIR}/Topology"
	"${PROJECT_SOURCE_DIR}/Geom"
	"${PROJECT_SOURCE_DIR}/SpecialMaps"
	"${PROJECT_SOURCE_DIR}/Refine"
	"${PROJECT_SOURCE_DIR}/Shape"
	"${PROJECT_SOURCE_DIR}/Material"
	"${PROJECT_SOURCE_DIR}/Multigrid"
	"${PROJECT_SOURCE_DIR}/Mesh"
	"${PROJECT_SOURCE_DIR}/StrMatrix"
	"${PROJECT_SOURCE_DIR}/Integral"
	"${PROJECT_SOURCE_DIR}/Frontal"
	"${PROJECT_SOURCE_DIR}/Pre"
	"${PROJECT_SOURCE_DIR}/Post"
	"${PROJECT_SOURCE_DIR}/Analysis"
	"${PROJECT_SOURCE_DIR}/LinearSolvers"
	"${PROJECT_SOURCE_DIR}/External"
	"${PROJECT_SOURCE_DIR}/External/sloan"
 )

configure_file(${PROJECT_SOURCE_DIR}/config.h.in ${PROJECT_BINARY_DIR}/config.h)
configure_file(${PROJECT_SOURCE_DIR}/PZ_USER.mk.in ${PROJECT_BINARY_DIR}/PZ_USER.mk)

add_subdirectory(Util)
add_subdirectory(Common)
add_subdirectory(Save)
add_subdirectory(Integral)
#add_subdirectory(LinearSolvers)
add_subdirectory(Matrix)
add_subdirectory(Topology)
add_subdirectory(Geom)
add_subdirectory(SpecialMaps)
add_subdirectory(Shape)
add_subdirectory(Refine)
add_subdirectory(External)
add_subdirectory(Material)
add_subdirectory(Mesh)
add_subdirectory(Analysis)
#add_subdirectory(Cosys)
add_subdirectory(Multigrid)
add_subdirectory(Post)
add_subdirectory(Frontal)
add_subdirectory(StrMatrix)
add_subdirectory(Pre)
add_subdirectory(SubStruct)
add_subdirectory(lib)
add_subdirectory(Projects)

enable_testing()

add_test (TestShapeRuns Projects/testshape/TestShape)

add_test (TestGeomRuns Projects/testgeom/TestGeom)

include(CTest)
########### install files ###############

#install(FILES  PZ_user.mk DESTINATION )



#original Makefile.am contents follow:

##
## Top Makefile.am for NeoPZ
##
## It will be processed by GNU automake to produce Makefile.in
##
## $Id: Makefile.am,v 1.10 2009-03-26 20:35:33 fortiago Exp $
##
#
#usermkdir = $(prefix)
#
#usermk_DATA = PZ_user.mk
#
#SUBDIRS = SpecialMaps Topology External Common Mesh Refine Geom Util Matrix Cosys Multigrid Material Post Shape Integral Analysis Frontal StrMatrix Pre LinearSolvers Save lib Projects
#
#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#PACKAGE = PZ
#DEVELOPER_ID = LabMeC
#
#clean-project:
#	cd Projects && make distclean
#
#backup: clean-project distclean
#	@ echo Making $(PACKAGE) backup...
#	@ echo -----------------------------------------------------------------
#	@ echo Generating tarball...
#	@ cd .. && tar --exclude={*.a,*.so,*.bz2,.tm_*,.snprj} -cf \
#          $(PACKAGE)-$(DEVELOPER_ID)-`date +%Y%m%d_%Hh%Mm`.tar ./Pantheon
#	@ echo Compacting...
#	@ cd .. && bzip2 $(PACKAGE)-*.tar
#	@ echo Generated file:
#	@ echo $(PACKAGE)-$(DEVELOPER_ID)-`date +%Y%m%d_%Hh%Mm`.tar
#	@ echo -----------------------------------------------------------------
#
### --| NeoPZ |------------------------------------------------------------------
#AM_CXXFLAGS = -D_AUTODIFF -DUSING_METIS
