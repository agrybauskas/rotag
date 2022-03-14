#ifndef _PDBXPARSER_H_
#define _PDBXPARSER_H_

#include <map>
#include <string>
#include <vector>

typedef std::map< std::string, std::map< std::string, std::vector<std::string> > > pdbx_loop_data;

std::vector<pdbx_loop_data> obtain_pdbx_loop( std::string pdbx_file,
                                              std::vector<std::string> categories,
                                              bool read_until_end );

#endif
