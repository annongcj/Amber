#CMake file which creates amber.sh and amber.csh, which are sourced by users to set all of the variables for using amber.
# Must be included after PythonConfig, PerlConfig, 3rdPartyTools, and Packaging.

message(STATUS "Generating amber source scripts")

# directory where all source scripts are stored
set(PATHSCRIPT_DIR ${CMAKE_CURRENT_LIST_DIR}/pathscripts)
set(PATHSCRIPT_OUTPUT_DIR ${CMAKE_BINARY_DIR}/pathscripts)

# variables for scripts to use (must convert from CMake logicals to 1 or 0)
# --------------------------------------------------------------------
if(BUILD_PERL)
	set(AMBERSH_PERL_SUPPORT 1)
else()
	set(AMBERSH_PERL_SUPPORT 0)
endif()

if(BUILD_PYTHON)
	set(AMBERSH_PYTHON_SUPPORT 1)
else()
	set(AMBERSH_PYTHON_SUPPORT 0)
endif()

# true if the platform supports LD_LIBRARY_PATH
if(TARGET_LINUX OR TARGET_OSX)
	set(AMBERSH_SUPPORTS_LDLIBPATH 1)
else()
	set(AMBERSH_SUPPORTS_LDLIBPATH 0)
endif()

if(BUILD_QUICK)
	set(AMBERSH_QUICK_SUPPORT 1)
else()
	set(AMBERSH_QUICK_SUPPORT 0)
endif()

string(TIMESTAMP CONFIGURATION_TIMESTAMP "%Y-%m-%d at %H:%M:%S")

# get_filename_component(MKL_LIBRARY_DIR ${MKL_CORE_LIBRARY} DIRECTORY)

# extra path and lib path folders
# --------------------------------------------------------------------

# --------------------------------------------------------------------
if(TARGET_WINDOWS)
	#create batch and powershell files
	
	# paths of files to be installed
	set(VAR_FILE_BAT ${PATHSCRIPT_OUTPUT_DIR}/amber.bat)
	set(VAR_FILE_PSH ${PATHSCRIPT_OUTPUT_DIR}/amber.ps1)
	set(INTERACTIVE_FILE_BAT ${PATHSCRIPT_DIR}/amber-interactive.bat)
	
	# Miniconda needs to be added to that path if we're using it
	set(EXTRA_PATH_FOLDERS "")
	
	if(DOWNLOAD_MINICONDA)
		set(EXTRA_PATH_FOLDERS "%AMBERHOME%\\miniconda;%AMBERHOME%\\miniconda\\Scripts;")
	endif()
	
	# MKL needs to be on the path so that it can dlopen itself
	# see cmake-buildscripts#95
	if(mkl_ENABLED)
		string(REPLACE "/" "\\" MKL_LIB_DIR_WIN "${MKL_LIBRARY_DIR}")
		set(EXTRA_PATH_FOLDERS "${EXTRA_PATH_FOLDERS}${MKL_LIB_DIR_WIN};")
	endif()
	
	# convert environment variable references (%FOO%) to Powershell syntax ($env:FOO)
	string(REGEX REPLACE "%([a-zA-Z]+)%" "$env:\\1" EXTRA_PATH_FOLDERS_PSH "${EXTRA_PATH_FOLDERS}")
	
	# perl path
	if(BUILD_PERL)
		# get windows path for perl lib
		string(REPLACE "/" "\\" PERL_MODULE_DIR_WIN ${PERL_MODULE_PATH})
	else()
		set(PERL_MODULE_DIR_WIN UNUSED)
	endif()

	# python path
	if(BUILD_PYTHON)
		string(REPLACE "/" "\\" RELATIVE_PYTHONPATH_WIN "${PREFIX_RELATIVE_PYTHONPATH}")
	else()
		set(RELATIVE_PYTHONPATH_WIN UNUSED)
	endif()
	
	configure_file(${PATHSCRIPT_DIR}/amber.in.bat ${VAR_FILE_BAT} @ONLY)
	configure_file(${PATHSCRIPT_DIR}/amber.in.ps1 ${VAR_FILE_PSH} @ONLY)
	
	#wrapper script which starts an interactive shell.
	install(PROGRAMS ${VAR_FILE_BAT} ${VAR_FILE_PSH} ${INTERACTIVE_FILE_BAT} DESTINATION "${CMAKE_INSTALL_POSTFIX}.")
endif()

# Generate bash scripts
# --------------------------------------------------------------------

# NOTE: Even on Windows, we still generate the Bash scripts because they might be useful to users using MSYS or Cygwin

# paths of files to be installed
set(VAR_FILE_SH ${PATHSCRIPT_OUTPUT_DIR}/amber.sh)
set(VAR_FILE_CSH ${PATHSCRIPT_OUTPUT_DIR}/amber.csh)
set(INTERACTIVE_FILE_SH ${PATHSCRIPT_DIR}/amber-interactive.sh)

if(TARGET_LINUX)
	set(LIB_PATH_VAR LD_LIBRARY_PATH)
	
	# determine folders to put on the library path
	set(LIB_PATH_DIRECTORIES "$AMBERHOME/lib")
	
	# MKL needs to be on the library path so that it can dlopen itself
	# see cmake-buildscripts#95
	if(mkl_ENABLED)
		set(LIB_PATH_DIRECTORIES "${LIB_PATH_DIRECTORIES}:${MKL_LIBRARY_DIR}")
	endif()
elseif(TARGET_OSX)
	set(LIB_PATH_VAR DYLD_FALLBACK_LIBRARY_PATH)

	# determine folders to put on the library path
	set(LIB_PATH_DIRECTORIES "$AMBERHOME/lib")
else()
	set(LIB_PATH_VAR UNUSED)
	set(LIB_PATH_DIRECTORIES UNUSED)
endif() 

set(EXTRA_PATH_PART "")
# if we're using minconda, we need to add that to the path if symlinks are not supported
if(DOWNLOAD_MINICONDA AND TARGET_WINDOWS)
	set(EXTRA_PATH_PART "$AMBERHOME/miniconda/bin:")
endif()


# perl path
if(BUILD_PERL)
	# get windows path for perl lib
	set(AMBERSH_PERL_MODULE_DIR ${PERL_MODULE_PATH})
else()
	set(AMBERSH_PERL_MODULE_DIR UNUSED)
endif()

# python path
if(BUILD_PYTHON)
	set(AMBERSH_RELATIVE_PYTHONPATH "${PREFIX_RELATIVE_PYTHONPATH}")
else()
	set(AMBERSH_RELATIVE_PYTHONPATH UNUSED)
endif()

configure_file(${PATHSCRIPT_DIR}/amber.in.sh ${VAR_FILE_SH} @ONLY)
configure_file(${PATHSCRIPT_DIR}/amber.in.csh ${VAR_FILE_CSH} @ONLY)
	
#put the scripts on the root dir of the install prefix
install(FILES  ${VAR_FILE_SH} ${VAR_FILE_CSH} DESTINATION "${CMAKE_INSTALL_POSTFIX}.")
install(PROGRAMS  ${INTERACTIVE_FILE_SH} DESTINATION "${CMAKE_INSTALL_POSTFIX}.")
