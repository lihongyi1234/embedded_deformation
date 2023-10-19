/*
 *   compute the normals using the faces
 *   by R. Falque
 *   16/01/2020
 */

#ifndef FOLDER_CHECKS_H
#define FOLDER_CHECKS_H

#include <dirent.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>
#include <string>

inline bool empty_folder(std::string folder_path) {
  // These are data types defined in the "dirent" header
  DIR *theFolder = opendir(folder_path.c_str());
  struct dirent *next_file;
  std::string filepath;

  while ((next_file = readdir(theFolder)) != NULL) {
    // build the path for each file in the folder
    filepath = folder_path + next_file->d_name;
    if (filepath.substr(filepath.size() - 4) == ".png")
      remove(filepath.c_str());
  }
  closedir(theFolder);

  return true;
};

inline bool does_folder_exist(std::string folder_path) {
  bool folder_exist;

  DIR *dir = opendir(folder_path.c_str());
  if (dir) {
    folder_exist = true;
    closedir(dir);
  } else if (ENOENT == errno) {
    folder_exist = false;
  } else {
    folder_exist = false;
  }

  return folder_exist;
};

inline bool create_folder(std::string folder_path) {
  mkdir(folder_path.c_str(), 0755);
  std::cout << "Progress: create the folder:" << folder_path << std::endl;
  return true;
}

#endif
