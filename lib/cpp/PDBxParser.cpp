#include <iostream>
#include <string>
#include <vector>
#include <boost/algorithm/string/join.hpp>
#include "PDBxParser.h"

/*
  Obtains pdbx loops for a specified categories.
  Input:
      pdbx_file - PDBx file path;
      categories - list of specified categories.
      read_until_end - boolean flag for reading whole pdbx file or stdin.
  Output:
      pdbx_loop_data - data structure for loop data or list of data structure.
*/

void obtain_pdbx_loop( std::string pdbx_file,
                       std::vector<std::string> categories,
                       bool read_until_end = false )
{
    std::vector<std::string> current_categories;
    std::vector<std::string> attributes;
    std::vector<std::string> data; /* Will be used for storing atom data temporarily.*/

    std::string category_regexp = boost::algorithm::join( categories , "|");
    bool is_reading_lines = false; /* Starts/stops reading lines at certain flags. */
}
