CMAKE_MINIMUM_REQUIRED ( VERSION 2.8 )
INCLUDE ( CheckFunctionExists )
INCLUDE ( CheckIncludeFiles ) 

# Disable the VC++ warnings that occur in the default configuration when compiling the fox toolkit.
IF ( MSVC )
  ADD_DEFINITIONS ( /wd4018 /wd4307 /wd4146 /wd4800 )
ENDIF ( MSVC )

# Silence all GCC warnings triggered by the fox toolkit.
IF ( CMAKE_COMPILER_IS_GNUCXX )
  ADD_DEFINITIONS ( -w )
  # For some reason march=native causes fox to crash when opening dialogs.
  STRING ( REPLACE "-march=native" "" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE} )
ENDIF ( CMAKE_COMPILER_IS_GNUCXX )

FILE ( GLOB HDRS include/*.h )
FILE ( GLOB SRCS src/*.cpp ) 

# Fox needs to be informed about the presence of certain headers and functions.
# The list below is most likely incomplete but at least sufficient for OS X 10.8 with MacPorts gcc48.
CHECK_INCLUDE_FILES ( dirent.h HAVE_DIRENT_H )
IF ( HAVE_DIRENT_H )
  ADD_DEFINITIONS ( -DHAVE_DIRENT_H=1 )
ENDIF ( )
CHECK_INCLUDE_FILES ( unistd.h HAVE_UNISTD_H )
IF ( HAVE_UNISTD_H )
  ADD_DEFINITIONS ( -DHAVE_UNISTD_H=1 )
ENDIF ( )
CHECK_FUNCTION_EXISTS ( strtoll HAVE_STRTOLL )
IF ( HAVE_STRTOLL )
  ADD_DEFINITIONS ( -DHAVE_STRTOLL=1 )
ENDIF ( )
CHECK_FUNCTION_EXISTS ( strtoull HAVE_STRTOULL )
IF ( HAVE_STRTOULL )
  ADD_DEFINITIONS ( -DHAVE_STRTOULL=1 )
ENDIF ( )
CHECK_FUNCTION_EXISTS ( vsscanf HAVE_VSSCANF )
IF ( HAVE_VSSCANF )
  ADD_DEFINITIONS ( -DHAVE_VSSCANF=1 )
  LIST ( REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/vsscanf.cpp )
ENDIF ( )
CHECK_INCLUDE_FILES ( time.h HAVE_TIME_H )
CHECK_INCLUDE_FILES ( sys/time.h HAVE_SYS_TIME_H )
IF ( HAVE_TIME_H AND HAVE_SYS_TIME_H )
  ADD_DEFINITIONS ( -DTIME_WITH_SYS_TIME=1 )
ENDIF ( )

# Some of the source files are not needed under some conditions.
IF ( HAVE_STRTOLL AND HAVE_STRTOULL )
  LIST ( REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/strtoll.cpp )
ENDIF ( )
IF ( UNIX )
  LIST ( REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/fxkeyboard.cpp )
ENDIF ( )

# The following two source files are only useful if fox is configured to use lbzip2 and zlib.
LIST ( REMOVE_ITEM SRCS ${CMAKE_CURRENT_SOURCE_DIR}/src/FXBZFileStream.cpp ${CMAKE_CURRENT_SOURCE_DIR}/src/FXGZFileStream.cpp )

INCLUDE_DIRECTORIES ( include )

SET ( EXTERNAL_FOX_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include CACHE INTERNAL "" )

ADD_LIBRARY ( FOX-1.6 STATIC ${SRCS} ${HDRS} )
IF ( APPLE )
  FIND_LIBRARY ( CORE_SERVICES_LIBRARY NAMES CoreServices )
  TARGET_LINK_LIBRARIES( FOX-1.6 ${CORE_SERVICES_LIBRARY} )
ENDIF ( )

SET ( EXTERNAL_FOX_LIBRARIES FOX-1.6 CACHE INTERNAL "" )