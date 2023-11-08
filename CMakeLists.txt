cmake_minimum_required(VERSION 3.2)
project(di2hap VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available

find_package(savvy REQUIRED)

add_executable(di2hap main.cpp)
target_compile_definitions(di2hap PUBLIC -DVERSION="${PROJECT_VERSION}")
target_link_libraries(di2hap savvy)

install(TARGETS di2hap RUNTIME DESTINATION bin)
