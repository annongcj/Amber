@echo off
rem Run this script to add the variables necessary to use Amber to your shell.
rem This script must be located in the Amber root folder!

rem Amber was configured on @CONFIGURATION_TIMESTAMP@
	
set AMBERHOME=%~dp0
set AMBERHOME=%AMBERHOME:~0,-1%

set "PATH=%AMBERHOME%\bin;@EXTRA_PATH_FOLDERS@%PATH%"

rem add Perl libs to path (only enabled if Amber was compiled with Perl support)
if @AMBERSH_PERL_SUPPORT@ equ 1 ( 
	if ["%PERL5LIB%"] equ [] (
		set "PERL5LIB=%AMBERHOME%\@PERL_MODULE_DIR_WIN@"
	) else (
		set "PERL5LIB=%AMBERHOME%\@PERL_MODULE_DIR_WIN@;%PERL5LIB%"
	)
)

rem add Python packages to path (only enabled if Amber was compiled with Python support)
if @AMBERSH_PYTHON_SUPPORT@ equ 1 ( 
	if ["%PYTHONPATH%"] equ [] (
		set "PYTHONPATH=%AMBERHOME%@RELATIVE_PYTHONPATH_WIN@"
	) else (
		set "PYTHONPATH=%AMBERHOME%@RELATIVE_PYTHONPATH_WIN@;%PYTHONPATH%"
	)
)

rem Tell QUICK where to find its data files
if @AMBERSH_QUICK_SUPPORT@ equ 1 ( 
	set "QUICK_BASIS=%AMBERHOME%\AmberTools\src\quick\basis"
)