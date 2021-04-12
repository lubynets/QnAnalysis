set(SOURCES Fitter.cpp
            ShapeContainer.cpp
            ShapeContainerTensor.cpp
            NDimFrame.cpp
            GraphExtractor.cpp
            ShapeFitter.cpp)
            
string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(QnAnalysisDiscriminator SHARED ${SOURCES} G__QnAnalysisDiscriminator.cxx)
ROOT_GENERATE_DICTIONARY(G__QnAnalysisDiscriminator ${HEADERS} LINKDEF QnAnalysisDiscriminatorLinkDef.h)
target_link_libraries(QnAnalysisDiscriminator ${ROOT_LIBRARIES} QnToolsBase)

add_subdirectory(exe)