#ifndef _PDBXPARSER_H_
#define _PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>

typedef std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::string>>>> pdbx_data;

void obtain_pdbx_data(std::string, std::vector<std::string>);

pdbx_data obtain_pdbx_line(std::string, std::vector<std::string>);

void obtain_pdbx_loop(std::string, std::vector<std::string>);

#endif
