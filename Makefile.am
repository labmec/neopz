#
# Top Makefile.am for NeoPZ
#
# It will be processed by GNU automake to produce Makefile.in
#
# $Id: Makefile.am,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $
#

usermkdir = $(prefix)

usermk_DATA = PZ_user.mk

SUBDIRS = Util Matrix Cosys Multigrid Material Post Shape Integral Geom \
	Mesh Refine External Analysis Frontal StrMatrix Pre Common \
	LinearSolvers lib

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
