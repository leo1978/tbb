cmake_minimum_required(VERSION 2.8)

find_package(TBB REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${TBB_INCLUDE_DIRS})
include_directories(${VTK_INCLUDE_DIRS})

add_executable(parallel_example parallel_example.cxx)
if (TARGET vtkHybrid)
  target_link_libraries(parallel_example ${TBB_LIBRARY} vtkHybrid)
else()
  target_link_libraries(parallel_example ${TBB_LIBRARY} vtkFiltersCore vtkImagingCore)
endif()

add_executable(atomic_example atomic_example.cxx)
if (TARGET vtkCommon)
  target_link_libraries(atomic_example ${TBB_LIBRARY} vtkCommon)
else()
  target_link_libraries(atomic_example ${TBB_LIBRARY} vtkCommonCore)
endif()

add_executable(register_thread_safety register_thread_safety.cxx)
if (TARGET vtkHybrid)
  target_link_libraries(register_thread_safety ${TBB_LIBRARY} vtkHybrid)
else()
  target_link_libraries(register_thread_safety ${TBB_LIBRARY} vtkFiltersCore vtkImagingCore)
endif()
