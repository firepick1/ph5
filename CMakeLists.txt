cmake_minimum_required(VERSION 2.8)

project( ph5 )

SET(PROJECT_VERSION_MAJOR 0)
SET(PROJECT_VERSION_MINOR 1)
set(PROJECT_VERSION_PATCH 0)
set(PROJECT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${PROJECT_VERSION_PATCH})


IF(WIN32)
  MESSAGE(STATUS "Detecting WINDOWS build")
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}")
  # By default, do not warn when built on machines using only VS Express:
  IF(NOT DEFINED CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS)
    SET(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS ON)
  ENDIF()
  SET(JANSSON_DIR "${CMAKE_SOURCE_DIR}/jansson/msvc" CACHE STRING "Jansson directory") 
  SET(JANSSON_INCLUDE_DIRS "${JANSSON_DIR}/include")
  SET(JANSSON_LIB "${JANSSON_DIR}/lib/Debug/jansson_d.lib")
  SET(JANSSON_DLL "${JANSSON_DIR}/bin/Debug/jansson_d.dll")
  MESSAGE(STATUS "JANSSON_DLL:${JANSSON_DLL}")
  MESSAGE(STATUS "JANSSON_LIB:${JANSSON_LIB}")
  SET(PH5_LIB "_ph5.lib")
  SET(PH5_DLL "_ph5.dll")
ELSE(WIN32)
  MESSAGE(STATUS "Detecting LINUX build")
  SET(JANSSON_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/jannson/src")
  SET(JANSSON_LIB "libjansson.so")
  SET(PH5_LIB "lib_ph5.so")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -g -Wno-format-extra-args")
  SET(CMAKE_SHARED_LINKER_FLAGS_DEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG} -g")
ENDIF(WIN32)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${CMAKE_SOURCE_DIR}/cmake/Modules/")

SET(COMPILE_DEFINITIONS -Werror)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_FILE_OFFSET_BITS=64")

get_property(dirs DIRECTORY . PROPERTY INCLUDE_DIRECTORIES)
message("INCLUDE_DIRECTORIES:${dirs}")

SET(BUILD_TARGET_DIR "${CMAKE_CURRENT_SOURCE_DIR}/target" CACHE STRING "Put all binary output in target")
SET(LIBRARY_OUTPUT_PATH "${BUILD_TARGET_DIR}")
SET(EXECUTABLE_OUTPUT_PATH "${BUILD_TARGET_DIR}")
SET(RUNTIME_OUTPUT_NAME ph5)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

include_directories( ${PROJECT_SOURCE_DIR} ${BUILD_TARGET_DIR} ${JANSSON_INCLUDE_DIRS} )
include_directories( /usr/local/include )
include_directories( /usr/include )

include (InstallRequiredSystemLibraries)
SET(CPACK_RESOURCE_FILE_LICENSE  
     "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_GENERATED_WARNING "GENERATED FILE!!! DO NOT CHANGE!!!")
SET(CPACK_PACKAGE_VERSION_MAJOR "${PROJECT_VERSION_MAJOR}")
SET(CPACK_PACKAGE_VERSION_MINOR "${PROJECT_VERSION_MINOR}")
include (CPack)

get_property(dirs DIRECTORY . PROPERTY INCLUDE_DIRECTORIES)
message("INCLUDE_DIRECTORIES:")
foreach(dir ${dirs})
  message(STATUS "  dir='${dir}'")
endforeach()

configure_file (
  "${PROJECT_SOURCE_DIR}/version.h.in"
  "${BUILD_TARGET_DIR}/version.h"
  )

message("CMAKE_PREFIX_PATH:${CMAKE_PREFIX_PATH}")

link_directories( ${BUILD_TARGET_DIR} /usr/lib /usr/local/lib )

set(PH5_LIB_FILES
  Complex.cpp
  PH5Curve.cpp
  FireLog.cpp
  )

add_library(_ph5 SHARED ${PH5_LIB_FILES})
target_link_libraries(_ph5 ${JANSSON_LIB} )
set_target_properties(_ph5 PROPERTIES 
    VERSION ${PROJECT_VERSION_STRING} 
    SOVERSION ${PROJECT_VERSION_MAJOR}
    )
if(WIN32)
  add_custom_command(TARGET _ph5 PRE_LINK
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        "${JANSSON_LIB}"      
        "${PROJECT_SOURCE_DIR}/target/Debug/" 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        "${JANSSON_DLL}"      
        "${PROJECT_SOURCE_DIR}/target/" 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        "${JANSSON_DLL}"      
        "${PROJECT_SOURCE_DIR}/target/Debug/" 
    )
  add_custom_command(TARGET _ph5 POST_BUILD    
    COMMAND cd ${PROJECT_SOURCE_DIR}/target/Debug
    COMMAND copy ${PH5_DLL} ..
    COMMAND copy ${PH5_DLL} ..
    COMMAND cd ${PROJECT_SOURCE_DIR}/msvc
    )
endif(WIN32)
  
add_executable(ph5 ph5.cpp)
add_dependencies(ph5 _ph5)
target_link_libraries(ph5 ${JANSSON_LIB} ${PH5_LIB} )
if(WIN32)
  add_custom_command(TARGET ph5 POST_BUILD    
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        "${PROJECT_SOURCE_DIR}/target/Debug/ph5.exe"      
        "${PROJECT_SOURCE_DIR}/target/ph5.exe" 
        )
endif(WIN32)

add_executable(test 
  test/test.cpp)

add_dependencies(test _ph5)
target_link_libraries(test ${JANSSON_LIB} ${PH5_LIB} )
if(WIN32)
  add_custom_command(TARGET test POST_BUILD    
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  
        "${PROJECT_SOURCE_DIR}/target/Debug/test.exe"      
        "${PROJECT_SOURCE_DIR}/target/test.exe" 
        )
endif(WIN32)
#
# Installation preparation.
#

# Allow the user to override installation directories.
set(PH5_INSTALL_LIB_DIR       lib CACHE PATH "Installation directory for libraries")
set(PH5_INSTALL_BIN_DIR       bin CACHE PATH "Installation directory for executables")
set(PH5_INSTALL_INCLUDE_DIR   include CACHE PATH "Installation directory for header files")


INSTALL(TARGETS _ph5 DESTINATION ${PH5_INSTALL_LIB_DIR})
INSTALL(TARGETS ph5 DESTINATION ${PH5_INSTALL_BIN_DIR})
INSTALL(FILES FireLog.h ph5.hpp DESTINATION ${PH5_INSTALL_INCLUDE_DIR})
