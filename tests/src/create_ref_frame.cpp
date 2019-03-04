#include "LinearAlgebraCpp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>

int main( int argc, char* argv[] ) {
    std::vector<double> mid_atom_coord(3);
    std::vector<double> up_atom_coord(3);
    std::vector<double> side_atom_coord(3);

    std::vector<double> coords_xyz;

    std::string line;
    std::ifstream fh( argv[1] );
    if( fh.is_open() ) {
        while( getline( fh, line ) ) {
            std::vector<std::string> line_split;
            boost::split( line_split, line, boost::is_any_of(" ")  );
            for( int i = 0; i < line_split.size(); i++ ) {
                if( line_split[i] != "" ) {
                    coords_xyz.push_back( std::stod( line_split[i] ) );
                }
            }
        }
        fh.close();
    }

    for( int i = 0; i < coords_xyz.size(); i++ ) {
        if( i / 3 == 0 ) { mid_atom_coord[i] = coords_xyz.at( i ); }
        if( i / 3 == 1 ) { up_atom_coord[i-3] = coords_xyz.at( i ); }
        if( i / 3 == 2 ) { side_atom_coord[i-6] = coords_xyz.at( i ); }
    }

    create_ref_frame( mid_atom_coord, up_atom_coord, side_atom_coord );
    return 0;
}
