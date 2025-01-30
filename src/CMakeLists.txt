FILE(GLOB SRCS *.cpp)
ADD_LIBRARY(Native SHARED ${SRCS})

TARGET_LINK_LIBRARIES(Native  PUBLIC ${LAPACKE_LIBRARIES} 
                             ${CBLAS_LIBRARIES} 
                             ${BLAS_LIBRARIES}
                             OpenMP::OpenMP_CXX)
# note: when the solver is used in CUPyDO, we must link explicitly 
#       to OpenMP so that CUPyDO is able to find OpenMP at runtime 
#       even if the fluid solver uses another threading library
