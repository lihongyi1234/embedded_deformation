/*
to add in CMakeList:

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.45.0 COMPONENTS filesystem)

include_directories(${filename} ${Boost_INCLUDE_DIRS} )
target_link_libraries(${filename} ${Boost_LIBRARIES} )

*/

#ifndef FOLDER_FOLDER_LIST_FILES_H
#define FOLDER_FOLDER_LIST_FILES_H

#include <string>
#include <vector>

int load_folder(std::string folder_path, std::vector<std::string>& filepaths,
                std::vector<std::string>& filenames);

int list_subfolders(std::string folder_path, std::vector<std::string>& subfolders);

#endif