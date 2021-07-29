set PATH=C:\NSRG\external\pc_win32\bin\
set TCL_LIBRARY=C:/nsrg/external/pc_win32/lib/tcl8.3/
set TK_LIBRARY=C:/nsrg/external/pc_win32/lib/tk8.3/
.\template_based_matching.exe %1
if not errorlevel 0 pause