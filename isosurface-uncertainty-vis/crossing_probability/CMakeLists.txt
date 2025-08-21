cmake_minimum_required(VERSION 3.23)

project(contour_variance CXX)

find_package(Viskores REQUIRED)
find_package(Boost 1.87 REQUIRED COMPONENTS random)
find_package (Eigen3 3.3 REQUIRED)

add_executable(crossing-probability compute-crossing-probability.cxx)

target_link_libraries(crossing-probability
  PRIVATE viskores::filter_contour viskores::io Boost::random Eigen3::Eigen
  )

viskores_add_target_information(crossing-probability
  DROP_UNUSED_SYMBOLS
  MODIFY_CUDA_FLAGS
  DEVICE_SOURCES MVGaussianTemp.hpp
  )
