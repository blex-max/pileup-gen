cmake_minimum_required (VERSION 3.22)

project(htsgen
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

option(MAKE_PY "Compile python bindings" OFF)
option(MAKE_EXE "Compile htsgen executable" OFF)
option(MAKE_TEST "Compile and run tests binary" OFF)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- find deps ---//
find_package(PkgConfig QUIET)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
    plog
    GIT_REPOSITORY https://github.com/SergiusTheBest/plog
    GIT_TAG        1.1.10
    GIT_SHALLOW    true
)
FetchContent_MakeAvailable(plog)

# ------- htslib -------
set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install. Must be >=1.14, but this is not enforced.")
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(_htslib_inc AND _htslib_lib)
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib>=1.14)
  if(HTSLIB_FOUND)
    message(STATUS "Located htslib>=1.14 install via pkg-config.")
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY to use user install."
  )
endif()
# //--- end deps ---//

# //--- testing ---//
if(MAKE_TEST)
  Include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  add_executable(test-htsgen
    tests/test.cpp
  )
  target_include_directories(test-htsgen PRIVATE src)
  target_link_libraries(test-htsgen PRIVATE
    Catch2::Catch2WithMain
    ${HTSLIB_TARGET}
  )
endif(MAKE_TEST)
# //--- end testing ---//

# compilation flags
add_library(common_compile_opts INTERFACE)
target_compile_options(common_compile_opts INTERFACE
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
    -Wsign-conversion
    -Wfloat-equal
    -Wnon-virtual-dtor
    -Wold-style-cast
    -Woverloaded-virtual
    -Wnull-dereference
    -Wdouble-promotion
    -Wformat=2
)

if(ENABLE_DEBUG_FLAGS)
  target_compile_options(common_compile_opts INTERFACE
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -O0 -g
  )
  target_link_options(common_compile_opts INTERFACE
    -fsanitize=address,undefined
  )
endif()


add_library(htsgen-core STATIC
  src/generate-pileup.cpp
)
target_include_directories(htsgen-core PUBLIC
  src
)
target_link_libraries(htsgen-core PUBLIC
  ${HTSLIB_TARGET}
)
target_link_libraries(htsgen-core PRIVATE
  plog::plog
)


if (MAKE_EXE)
  FetchContent_Declare(
      argparse
      GIT_REPOSITORY https://github.com/p-ranav/argparse.git
      GIT_TAG v3.2
  )
  FetchContent_MakeAvailable(argparse)
  add_executable(htsgen
    src/main.cpp
  )
  target_link_libraries(htsgen PRIVATE
    common_compile_opts
    htsgen-core
    argparse
    plog::plog
  )
endif()


# link bindings
if(MAKE_PY)
  set(PYBIND11_FINDPYTHON NEW)
  FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v3.0.3
  )
  FetchContent_MakeAvailable(pybind11)
  pybind11_add_module(pyhtsgen
    python/bindings.cpp
  )
  target_link_libraries(pyhtsgen PRIVATE
    common_compile_opts
    htsgen-core
  )
  set_target_properties(pyhtsgen PROPERTIES
    OUTPUT_NAME htsgen
  )
  # for pip install
  install(TARGETS pyhtsgen
    LIBRARY DESTINATION .
  )
endif()


