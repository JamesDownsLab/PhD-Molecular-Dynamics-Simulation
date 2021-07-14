cmake_minimum_required(VERSION 3.19)
project(3DMolecularDynamics)

set(CMAKE_CXX_STANDARD 20)

add_executable(3DMolecularDynamics main.cpp Engine.cpp Engine.h Particle.cpp Particle.h Vector.cpp Vector.h setup.h)