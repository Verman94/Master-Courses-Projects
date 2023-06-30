# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=A_ADI - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to A_ADI - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "A_ADI - Win32 Release" && "$(CFG)" != "A_ADI - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "A_ADI.mak" CFG="A_ADI - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "A_ADI - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "A_ADI - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "A_ADI - Win32 Debug"
RSC=rc.exe
F90=fl32.exe

!IF  "$(CFG)" == "A_ADI - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\A_ADI.exe"

CLEAN : 
	-@erase ".\A_ADI.exe"
	-@erase ".\CAVITY.obj"
	-@erase ".\A_ADI.obj"
	-@erase ".\grid.obj"
	-@erase ".\pcol1.obj"

# ADD BASE F90 /Ox /c /nologo
# ADD F90 /Ox /c /nologo
F90_PROJ=/Ox /c /nologo 
# ADD BASE RSC /l 0x818 /d "NDEBUG"
# ADD RSC /l 0x818 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/A_ADI.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/A_ADI.pdb" /machine:I386 /out:"$(OUTDIR)/A_ADI.exe" 
LINK32_OBJS= \
	"$(INTDIR)/CAVITY.obj" \
	"$(INTDIR)/A_ADI.obj" \
	"$(INTDIR)/grid.obj" \
	"$(INTDIR)/pcol1.obj"

"$(OUTDIR)\A_ADI.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "A_ADI - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
OUTDIR=.
INTDIR=.

ALL : "$(OUTDIR)\A_ADI.exe"

CLEAN : 
	-@erase ".\A_ADI.exe"
	-@erase ".\CAVITY.obj"
	-@erase ".\A_ADI.obj"
	-@erase ".\grid.obj"
	-@erase ".\pcol1.obj"
	-@erase ".\A_ADI.ilk"
	-@erase ".\A_ADI.pdb"

# ADD BASE F90 /Zi /c /nologo
# ADD F90 /Zi /c /nologo
F90_PROJ=/Zi /c /nologo /Fd"A_ADI.pdb" 
# ADD BASE RSC /l 0x818 /d "_DEBUG"
# ADD RSC /l 0x818 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/A_ADI.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/A_ADI.pdb" /debug /machine:I386 /out:"$(OUTDIR)/A_ADI.exe" 
LINK32_OBJS= \
	"$(INTDIR)/CAVITY.obj" \
	"$(INTDIR)/A_ADI.obj" \
	"$(INTDIR)/grid.obj" \
	"$(INTDIR)/pcol1.obj"

"$(OUTDIR)\A_ADI.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for.obj:
   $(F90) $(F90_PROJ) $<  

.f.obj:
   $(F90) $(F90_PROJ) $<  

.f90.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "A_ADI - Win32 Release"
# Name "A_ADI - Win32 Debug"

!IF  "$(CFG)" == "A_ADI - Win32 Release"

!ELSEIF  "$(CFG)" == "A_ADI - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\A_ADI.f90

"$(INTDIR)\A_ADI.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE="C:\Users\Mohammad\Desktop\New folder\A_ADI.f90"

"$(INTDIR)\A_ADI.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE="\Desc\Book\Semester 2\CFD\Project#4\CODE\CAVITY\CAVITY.f90"

"$(INTDIR)\CAVITY.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE="C:\Users\Mohammad\Desktop\Laminar flow\grid.f"

"$(INTDIR)\grid.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
################################################################################
# Begin Source File

SOURCE="C:\Users\Mohammad\Desktop\Laminar flow\pcol1.f"

"$(INTDIR)\pcol1.obj" : $(SOURCE) "$(INTDIR)"
   $(F90) $(F90_PROJ) $(SOURCE)


# End Source File
# End Target
# End Project
################################################################################
