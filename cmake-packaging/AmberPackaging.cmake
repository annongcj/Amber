# Top-level package script that calls all the others.

if(AMBERTOOLS_ONLY)
	set(PACKAGE_NAME "Amber Tools")
	set(PACKAGE_FILENAME "AmberTools")
else()
	set(PACKAGE_NAME "Amber MD")
	set(PACKAGE_FILENAME "Amber")
endif()
set(ICO_ICON ${CMAKE_CURRENT_LIST_DIR}/amber_logo_no_background.ico)
set(ICO_UNINSTALL_ICON ${CMAKE_CURRENT_LIST_DIR}/amber_uninstall_no_background.ico)
set(ICNS_ICON ${CMAKE_CURRENT_LIST_DIR}/amber_logo_no_background.icns)

if(TARGET_OSX)
	set(STARTUP_FILE ${CMAKE_CURRENT_LIST_DIR}/osx-startup-script.sh)
elseif(TARGET_WINDOWS)
	set(STARTUP_FILE ${CMAKE_CURRENT_LIST_DIR}/pathscripts/amber-interactive.bat)
endif()

set(BUNDLE_IDENTIFIER org.ambermd.amber)
set(BUNDLE_SIGNATURE AMBR)
include(${CMAKE_CURRENT_LIST_DIR}/AmberComponents.cmake)
include(Packaging)