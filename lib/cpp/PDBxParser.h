#ifndef _PDBXPARSER_H_
#define _PDBXPARSER_H_

#include <string>
#include <vector>

void obtain_pdbx_loop( std::string pdbx_file,
                       std::vector<std::string> categories,
                       bool read_until_end );

#endif
