# Microsoft Developer Studio Project File - Name="pzvc" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=pzvc - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "pzvc.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "pzvc.mak" CFG="pzvc - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "pzvc - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "pzvc - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "pzvc - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "pzvc___Win32_Release"
# PROP BASE Intermediate_Dir "pzvc___Win32_Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /D "PZENVIRONMENT" /YX /FD /TP /c
# ADD BASE RSC /l 0x416 /d "NDEBUG"
# ADD RSC /l 0x416 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386

!ELSEIF  "$(CFG)" == "pzvc - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "pzvc___Win32_Debug"
# PROP BASE Intermediate_Dir "pzvc___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /D "PZENVIRONMENT" /FR /YX /FD /GZ /TP /c
# ADD BASE RSC /l 0x416 /d "_DEBUG"
# ADD RSC /l 0x416 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept

!ENDIF 

# Begin Target

# Name "pzvc - Win32 Release"
# Name "pzvc - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "material"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\material\pzelasmat.c
# End Source File
# Begin Source File

SOURCE=..\..\material\pzmat1dlin.c
# End Source File
# Begin Source File

SOURCE=..\..\material\pzmat2dlin.c
# End Source File
# Begin Source File

SOURCE=..\..\material\pzmaterial.c
# End Source File
# Begin Source File

SOURCE=..\..\material\pzmatplaca.c
# End Source File
# Begin Source File

SOURCE=..\..\material\pzmattest.c
# End Source File
# End Group
# Begin Group "matrix"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\matrix\pzblock.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzbndmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzespmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzfmatrix.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzlink.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzmatred.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzmatrix.c

!IF  "$(CFG)" == "pzvc - Win32 Release"

# SUBTRACT CPP /O<none>

!ELSEIF  "$(CFG)" == "pzvc - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzsbndmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzsespmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzsfulmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzshtmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzskylmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzsolve.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzstencil.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzsysmp.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pztempmat.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzworkpool.c
# End Source File
# Begin Source File

SOURCE=..\..\matrix\pzysmp.c
# End Source File
# End Group
# Begin Group "util"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\util\pzadmchunk.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzavlmap.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzavlmap_old.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzchunk.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzmanvector.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzmap.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzstack.c
# End Source File
# Begin Source File

SOURCE=..\..\util\pzvec.c
# End Source File
# End Group
# Begin Group "mesh"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\mesh\Elc1dgd.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzcmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzcompel.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzconnect.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzcosys.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelc1d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelcpoint.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelcq2d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelct2d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelg1d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelgpoint.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelgq2d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelgt2d.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzelmat.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzgeoel.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzgmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzgnode.c
# End Source File
# Begin Source File

SOURCE=..\..\mesh\pzintel.c
# End Source File
# End Group
# Begin Group "post"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\post\pzdxmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphel.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphel1d.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphel1dd.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphelq2d.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphelq2dd.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzgraphnode.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzmvmesh.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pztrigraph.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pztrigraphd.c
# End Source File
# Begin Source File

SOURCE=..\..\post\pzv3dmesh.c
# End Source File
# End Group
# Begin Group "multigrid"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\multigrid\pztransfer.c
# End Source File
# Begin Source File

SOURCE=..\..\multigrid\pztrnsform.c
# End Source File
# End Group
# Begin Group "analysis"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\analysis\pzanalysis.c
# End Source File
# Begin Source File

SOURCE=..\..\analysis\pzanalysiserror.c
# End Source File
# End Group
# Begin Group "integral"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\integral\pzquad.c
# End Source File
# End Group
# Begin Group "metis"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\..\external\metis\balance.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\bucketsort.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\ccgraph.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\coarsen.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\compress.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\debug.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\estmem.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\fm.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\fortran.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\frename.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\graph.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\initpart.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kvmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kwayfm.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kwayrefine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kwayvolfm.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\kwayvolrefine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\match.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mbalance.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mbalance2.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mcoarsen.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\memory.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mesh.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\meshpart.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mfm.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mfm2.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mincover.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\minitpart.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\minitpart2.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mkmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mkwayfmh.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mkwayrefine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mmatch.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mmd.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mpmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mrefine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mrefine2.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\mutil.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\myqsort.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\ometis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\parmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\pmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\pqueue.c
# End Source File
# Begin Source File

SOURCE=..\..\external\pzmetis.c
# End Source File
# Begin Source File

SOURCE=..\..\external\pzrenumbering.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\refine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\separator.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\sfm.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\srefine.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\stat.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\subdomains.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\timing.c
# End Source File
# Begin Source File

SOURCE=..\..\external\metis\util.c
# End Source File
# End Group
# Begin Source File

SOURCE=.\chapeu.c

!IF  "$(CFG)" == "pzvc - Win32 Release"

# PROP Exclude_From_Build 1

!ELSEIF  "$(CFG)" == "pzvc - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\Flavio.c

!IF  "$(CFG)" == "pzvc - Win32 Release"

# SUBTRACT CPP /WX

!ELSEIF  "$(CFG)" == "pzvc - Win32 Debug"

# PROP Exclude_From_Build 1

!ENDIF 

# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# End Group
# Begin Group "Resource Files"

# PROP Default_Filter "ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe"
# End Group
# End Target
# End Project
