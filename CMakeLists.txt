cmake_minimum_required(VERSION 3.19)
project(3DMolecularDynamics)

set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 3.3 REQUIRED)

add_executable(3DMolecularDynamics main.cpp Engine.cpp Engine.h Particle.cpp Particle.h BasePlate.cpp BasePlate.h Options.h Options.cpp nanoflann.h KDTreeVectorOfVectorsAdaptor.h)

target_link_libraries(3DMolecularDynamics Eigen3::Eigen)