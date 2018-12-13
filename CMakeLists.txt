cmake_minimum_required(VERSION 2.8.3)
project(MarkerBA)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
SET(CMAKE_BUILD_TYPE Debug)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/CMakeModules/")

find_package(OpenCV REQUIRED)
find_package(Eigen REQUIRED)
FIND_PACKAGE(Sophus REQUIRED)

FIND_PACKAGE(G2O REQUIRED)
INCLUDE_DIRECTORIES(
  ${G2O_INCLUDE_DIR}
  /usr/include/suitesparse # for cholmod
)
LIST(APPEND G2O_LIBS
  g2o_core  g2o_solver_csparse  g2o_csparse_extension g2o_types_sba g2o_solver_dense g2o_stuff g2o_parser g2o_solver_pcg  cholmod cxsparse
)

INCLUDE_DIRECTORIES(
  ${PROJECT_SOURCE_DIR}/include
  ${Eigen_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${Sophus_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} SHARED
  src/camera_model.cpp
  src/g2otypes_marker.cpp
)

#message(${PROJECT_SOURCE_DIR})
target_link_libraries(${PROJECT_NAME}
  ${OpenCV_LIBS}
  ${EIGEN3_LIBS}
  ${Sophus_LIBRARIES}
  ${Boost_LIBRARIES}
  ${G2O_LIBS}

)

add_executable(test_g2o_marker test/test_g2o_marker.cpp)
target_link_libraries(test_g2o_marker ${PROJECT_NAME})


