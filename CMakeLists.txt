cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_MODULE_PATH})
message(${CMAKE_MODULE_PATH})

find_package(TBB REQUIRED)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

include_directories(${TBB_INCLUDE_PATH})
include_directories(${VTK_INCLUDE_DIRS})

add_executable(parallel_example parallel_example.cxx)
target_link_libraries(parallel_example ${TBB_LIBRARY} vtkFiltersCore vtkImagingCore)

add_executable(parallel_example2 parallel_example2.cxx)
target_link_libraries(parallel_example2 ${TBB_LIBRARY} vtkFiltersCore vtkImagingCore)

add_executable(contouramr contouramr.cxx)
target_link_libraries(contouramr ${TBB_LIBRARY}  vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(serial-contouramr serial-contouramr.cxx)
target_link_libraries(serial-contouramr ${TBB_LIBRARY}  vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(threaded_streamtrace threaded_streamtrace.cxx)
target_link_libraries(threaded_streamtrace ${TBB_LIBRARY} vtkIOLegacy vtkImagingGeneral vtkFiltersFlowPaths vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(pcontour pcontour.cxx)
target_link_libraries(pcontour ${TBB_LIBRARY} vtkFiltersGeneral vtkIOLegacy vtkImagingGeneral vtkFiltersFlowPaths vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(arraytest arraytest.cxx)
target_link_libraries(arraytest ${TBB_LIBRARY} vtkFiltersGeneral vtkIOLegacy vtkImagingGeneral vtkFiltersFlowPaths vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(atomic_example atomic_example.cxx)
target_link_libraries(atomic_example ${TBB_LIBRARY} vtkCommonCore)

add_executable(register_thread_safety register_thread_safety.cxx)
target_link_libraries(register_thread_safety ${TBB_LIBRARY} vtkFiltersCore vtkImagingCore)

add_executable(parallel_amr parallel_amr.cxx)
target_link_libraries(parallel_amr ${TBB_LIBRARY}  vtkIOAMR vtkIOXML vtkFiltersAMR vtkRenderingVolumeOpenGL vtkRenderingFreeTypeOpenGL vtkInteractionStyle)

add_executable(testparallelfor testparallelfor.cxx)
target_link_libraries(testparallelfor ${TBB_LIBRARY} vtkFiltersCore pthread)

