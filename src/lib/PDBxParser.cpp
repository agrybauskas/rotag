#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/regex.hpp>
#include "PDBxParser.h"

/*
  Obtains pdbx data for the specified categories or items.
  Input:
      pdbx_file - PDBx file path;
      data_identifier - list of categories or items;
      read_stream - reads stream of pdbx files.
  Output:
      pdbx_data - data structure for pdbx data.
*/

void obtain_pdbx_data(std::string pdbx_file,
                      std::vector<std::string> data_identifier) {
  std::vector<std::string> pdbx_data;

  if(data_identifier.size() > 0) {
    std::vector<std::string> pdbxs;

    // Slurp whole pdbx file/files.
    // TODO: check if slurps multiple PDBxs in stream.
    int pdbx_counter = 0;
    std::string pdbx_line;
    std::ifstream fh(pdbx_file);
    if(fh.is_open()) {
      while(getline(fh, pdbx_line, '\0')){
        std::vector<std::string> pdbx_split;
        boost::algorithm::split_regex(pdbx_split, pdbx_line,
                                      boost::regex("(?<!\\S)data_"));
        for(int i = 0; i < pdbx_split.size(); i++) {
          if(pdbx_split[i] != "") {
              pdbxs.push_back("data_" + pdbx_split[i]);
          }
        }
      }
    }
    fh.close();

    // TODO: add warning.
    if(pdbx_counter == 0) {
    }

    for(int i; i < pdbxs.size(); i++) {
      obtain_pdbx_line(pdbxs[i], data_identifier);
    }
  }
}

/*
  Obtains pdbx lines for a specified items.
  Input:
      pdbx_file - PDBx file path;
      items - list of specified items.
  Output:
      pdbx_line_data - data structure for list of pdbx data structure.
*/

void obtain_pdbx_line(std::string pdbx_file,
                      std::vector<std::string> items) {
    std::string item_regexp = boost::algorithm::join(items, "|");
    boost::regex single_line_re{"(" + item_regexp + "|" + item_regexp +
                                ".\\S+)\\s+(?!;)('.+'|\\S+)"};
    boost::regex multi_line_re{"(" + item_regexp + "|" + item_regexp +
                               ".\\S+)\\s+(\n;[^;]+;)"};
    boost::smatch matches;
    while(boost::regex_search(pdbx_file, matches, single_line_re)) {
      std::cout << matches.str(1) << std::endl;
      pdbx_file = matches.suffix().str();
    }
}

// /*
//   Obtains pdbx loops for a specified categories.
//   Input:
//       pdbx_file - PDBx file path;
//       categories - list of specified categories.
//   Output:
//       pdbx_loop_data - data structure for loop data or list of data structure.
// */

// std::vector<pdbx_loop_data> obtain_pdbx_loop(std::string pdbx_file,
//                                              std::vector<std::string> categories)
// {
//     // NOTE: do not forget "ignore_missing_categories" option.

//     std::vector<std::string> category_list;
//     std::vector<std::string> attributes;
//     std::vector<std::string> data;

//     boost::regex category_regex(boost::algorithm::join(categories , "|"));

//     std::cout << category_regex << std::endl;
//     //NOTE: not sure if s/\[/\\[/g and s/\]/\\]/g is needed.
//     bool is_reading_lines = false;// Starts/stops reading lines at certain flags.

//     // std::vector<std::vector<std::string>> current_categories;
//     // std::vector<std::vector< std::vector<std::string>>> attributes;
//     // /* Will be used for storing atom data temporarily.*/
//     // std::vector<std::vector<std::vector<std::string>>> data;

//     // boost::regex category_regex(
//     //     "(" + boost::algorithm::join(categories , "|") + ")[.](.+)\n?$"
//     // );
//     // boost::regex data_regex("^data_");

//     // boost::regex reading_end_regex("^_|loop_|#");

//     // bool is_reading_lines = false; /* Starts/stops reading lines at certain flags. */

//     // std::ifstream fh(pdbx_file);
//     // std::string line;
//     // boost::smatch match;
//     // if(fh.is_open()) {
//     //     while(getline(fh, line)) {
//     //         if(boost::regex_search(line, match, data_regex) ||
//     //            current_categories.empty()) {
//     //             current_categories.push_back(std::vector<std::string>());
//     //             attributes.push_back(std::vector<std::vector<std::string>>());
//     //             data.push_back(std::vector<std::vector<std::string>>());
//     //         } else if(boost::regex_search(line, match, category_regex)) {
//     //             if(current_categories.back().empty() ||
//     //                current_categories.back().back() != match[1]) {
//     //                 std::string category(match[1]);
//     //                 current_categories.back().push_back(category);
//     //                 attributes.back().push_back(std::vector<std::string>());
//     //                 data.back().push_back(std::vector<std::string>());
//     //             }
//     //             std::string attribute(match[2]);
//     //             boost::trim(attribute);
//     //             attributes.back().back().push_back(attribute);
//     //             is_reading_lines = true;
//     //         } else if(is_reading_lines &&
//     //                   boost::regex_search(line, match, reading_end_regex)) {
//     //             if(current_categories.size() == current_categories.size() &&
//     //                ! read_until_end) {
//     //                 break;
//     //             }
//     //             is_reading_lines = false;
//     //         } else if(is_reading_lines) {
//     //             std::vector<std::string> data_split;
//     //             boost::split(data_split, line, boost::is_any_of(" "));
//     //             for(int i = 0; i < data_split.size(); i++) {
//     //                 if(data_split[i] != "") {
//     //                     data.back().back().push_back(data_split[i]);
//     //                 }
//     //             }
//     //         }
//     //     }
//     //     fh.close();
//     // }

//     // /* Generates hash from three lists. */
//     std::vector<pdbx_loop_data> pdbx_loop_data_list;
//     // for(int i = 0; i < current_categories.size(); i++) {
//     //     pdbx_loop_data pdbx_loop_data;
//     //     for(int j = 0; j < current_categories.back().size(); j++) {
//     //         pdbx_loop_data[current_categories[i][j]]["attributes"] = attributes[i][j];
//     //         pdbx_loop_data[current_categories[i][j]]["data"] = data[i][j];
//     //     }

//     //     pdbx_loop_data_list.push_back(pdbx_loop_data);
//     // }

//     return pdbx_loop_data_list;
// }
