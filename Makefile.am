#
# Top Makefile.am for NeoPZ
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id: Makefile.am,v 1.10 2009-03-26 20:35:33 fortiago Exp $
#

usermkdir = $(prefix)

usermk_DATA = PZ_user.mk

SUBDIRS = SpecialMaps Topology External Common Mesh Refine Geom Util Matrix Cosys Multigrid Material Post Shape Integral Analysis Frontal StrMatrix Pre LinearSolvers Save lib Projects

DISTCLEANFILES = @PZ_DISTCLEANFILES@
CLEANFILES     = @PZ_CLEANFILES@

PACKAGE = PZ
DEVELOPER_ID = LabMeC

clean-project:
	cd Projects && make distclean

backup: clean-project distclean
	@ echo Making $(PACKAGE) backup...
	@ echo -----------------------------------------------------------------
	@ echo Generating tarball...
	@ cd .. && tar --exclude={*.a,*.so,*.bz2,.tm_*,.snprj} -cf \
          $(PACKAGE)-$(DEVELOPER_ID)-`date +%Y%m%d_%Hh%Mm`.tar ./Pantheon
	@ echo Compacting...
	@ cd .. && bzip2 $(PACKAGE)-*.tar
	@ echo Generated file:
	@ echo $(PACKAGE)-$(DEVELOPER_ID)-`date +%Y%m%d_%Hh%Mm`.tar
	@ echo -----------------------------------------------------------------

## --| NeoPZ |------------------------------------------------------------------
AM_CXXFLAGS = -D_AUTODIFF -DUSING_METIS
