cmake_minimum_required(VERSION 3.8)
project(particle_sim CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED)
find_package(Ascent REQUIRED)

include_directories(${MPI_CXX_INCLUDE_PATH}
                    ${CONDUIT_INCLUDE_DIRS}
                    ${ASCENT_INCLUDE_DIRS})

add_executable(particle_sim particle_sim.cpp)
target_link_libraries(particle_sim ${MPI_CXX_LIBRARIES} ascent::ascent_mpi)

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/ascent_actions.yaml"
    "${CMAKE_CURRENT_BINARY_DIR}/ascent_actions.yaml"
    COPYONLY
)

