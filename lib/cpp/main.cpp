#include "PDBxParser.h"
#include "vector"

int main()
{
    std::vector<std::string> categories;
    categories.push_back( "_atom_site" );
    std::vector<std::string> unique_keys;
    unique_keys.push_back( "id" );
    pdbx_loop_unique( obtain_pdbx_loop( "serine-001.cif", categories, false ),
                      unique_keys, false );
    return 0;
}
