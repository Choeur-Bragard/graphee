cmake_minimum_required(VERSION 3.10)
project(graphee)

add_subdirectory(vendor/onelogger)

enable_testing()

set(CMAKE_CXX_STANDARD 14)
set(THREADS_PREFER_PTHREAD_FLAG ON)

include_directories(include)
include_directories(vendor/onelogger/onemonitor/vendor/prometheus-cpp/3rdparty/civetweb/include)
include_directories(vendor/onelogger/onemonitor/vendor/prometheus-cpp/core/include)
include_directories(vendor/onelogger/onemonitor/vendor/prometheus-cpp/pull/include)
include_directories(vendor/onelogger/include)
include_directories(vendor/onelogger/onemonitor/include)
include_directories(src/snappy/build/)


link_directories(/usr/lib
  /usr/local/lib
  build/vendor/onelogger/onemonitor/vendor/prometheus-cpp/pull
  build/vendor/onelogger/onemonitor/vendor/prometheus-cpp/core
  src/snappy/build/ 
  )

find_package(Threads REQUIRED)
find_package(Boost COMPONENTS unit_test_framework REQUIRED)

set(SOURCES
  src/pagerank.cpp
  src/edgelist.cpp
  src/utils.cpp
)

set(LIBS
  prometheus-cpp-pull
  onemonitor
  onelogger
  snappy
  omp
)
  

if(APPLE)
  set(LIBS ${LIBS} boost_thread-mt)
else()
  set(LIBS ${LIBS} boost_thread)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g ")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall -O3 ")

add_executable(graphee ${SOURCES})
target_link_libraries(graphee ${LIBS})

# declares a test with our executable
add_test(NAME tests COMMAND tests)

install(TARGETS graphee
        DESTINATION /build/)

