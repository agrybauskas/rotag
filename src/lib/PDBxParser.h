#ifndef _PDBXPARSER_H_
#define _PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>

// typedef std::map<std::string, std::map<std::string, std::vector<std::string>>> pdbx_data;

// std::vector<pdbx_loop_data> obtain_pdbx_loop(std::string pdbx_file,
//                                              std::vector<std::string> categories);

void obtain_pdbx_data(std::string pdbx_file,
                      std::vector<std::string> data_identifier,
                      bool read_stream);

#endif
