set(BUILD_SAMPLES OFF)
set(BUILD_EXAMPLES OFF)
set(BUILD_DEPS ON)
set(BUILD_SHARED_LIBS OFF)
add_subdirectory(or-tools)

set(CXXOPTS_BUILD_EXAMPLES OFF)
set(CXXOPTS_BUILD_TESTS OFF)
set(CXXOPTS_ENABLE_INSTALL OFF)
set(CXXOPTS_USE_UNICODE_HELP OFF)
add_subdirectory(cxxopts)