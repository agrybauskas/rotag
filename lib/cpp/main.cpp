#include "PDBxParser.h"
#include "vector"

int main()
{
    std::vector<std::string> categories;
    categories.push_back( "_atom_site" );
    obtain_pdbx_loop( "serine-001.cif", categories, false );
    return 0;
}
