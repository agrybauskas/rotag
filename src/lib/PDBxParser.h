#ifndef _PDBXPARSER_H_
#define _PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>
#include <boost/any.hpp>

struct pdbx_category {
  std::vector<std::string>  attributes = {};
  bool is_loop = false;
  bool is_indexed = false;

  boost::any data = {};
};

typedef std::map<std::string, pdbx_category> pdbx_data;

void obtain_pdbx_data(std::string, std::vector<std::string>);

pdbx_data obtain_pdbx_line(std::string, std::vector<std::string>);

void obtain_pdbx_loop(std::string, std::vector<std::string>);
#endif
