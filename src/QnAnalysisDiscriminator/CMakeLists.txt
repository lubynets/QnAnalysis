# include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_INCLUDE_DIRECTORIES})

add_executable(QnAnalysisDiscriminator main.cxx)
target_link_libraries(QnAnalysisDiscriminator ${ROOT_LIBRARIES})