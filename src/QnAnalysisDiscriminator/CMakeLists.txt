set(SOURCES Fitter.cpp
            ShapeContainer.cpp
            ShapeContainerTensor.cpp
            NDimFrame.cpp
            GraphExtractor.cpp)
            
string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(QnAnalysisDiscriminator SHARED ${SOURCES} G__QnAnalysisDiscriminator.cxx)
ROOT_GENERATE_DICTIONARY(G__QnAnalysisDiscriminator ${HEADERS} LINKDEF QnAnalysisDiscriminatorLinkDef.h)
target_link_libraries(QnAnalysisDiscriminator ${ROOT_LIBRARIES} QnToolsBase)

add_executable(discriminator discriminator.cxx)
add_dependencies(discriminator QnAnalysisDiscriminator)
target_link_libraries(discriminator QnAnalysisDiscriminator)

add_executable(shapewriter shapewriter.cxx)
add_dependencies(shapewriter QnAnalysisDiscriminator)
target_link_libraries(shapewriter QnAnalysisDiscriminator)

add_executable(shapetensorwriter shapetensorwriter.cxx)
add_dependencies(shapetensorwriter QnAnalysisDiscriminator)
target_link_libraries(shapetensorwriter QnAnalysisDiscriminator)

add_executable(discriminatortensor discriminatortensor.cxx)
add_dependencies(discriminatortensor QnAnalysisDiscriminator)
target_link_libraries(discriminatortensor QnAnalysisDiscriminator)

add_executable(mcfitter mcfitter.cxx)
add_dependencies(mcfitter QnAnalysisDiscriminator)
target_link_libraries(mcfitter QnAnalysisDiscriminator)

add_executable(mcextractor mcextractor.cxx)
add_dependencies(mcextractor QnAnalysisDiscriminator)
target_link_libraries(mcextractor QnAnalysisDiscriminator)
