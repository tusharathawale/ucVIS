cmake_minimum_required(VERSION 3.23)

project(contour_variance CXX)

find_package(Viskores REQUIRED)
find_package(Boost 1.87 REQUIRED COMPONENTS random)
find_package (Eigen3 3.3 REQUIRED)

add_executable(run-contour run-contour.cxx ContourVariance.cxx uniform_even_z_density.cxx uniform_kernel_polynomial.cxx)

target_link_libraries(run-contour
  PRIVATE viskores::filter_contour viskores::io Boost::random Eigen3::Eigen
  )

viskores_add_target_information(run-contour
  DROP_UNUSED_SYMBOLS
  MODIFY_CUDA_FLAGS
  DEVICE_SOURCES ContourVariance.cxx
  )
