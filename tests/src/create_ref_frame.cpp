#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include "LinearAlgebraCpp.h"

int main( int argc, char* argv[] ) {
    std::vector<double> mid_atom_coord;
    std::vector<double> up_atom_coord;
    std::vector<double> side_atom_coord;

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
        if( i / 3 == 0 ) { mid_atom_coord.push_back( coords_xyz.at( i ) ); }
        if( i / 3 == 1 ) { up_atom_coord.push_back( coords_xyz.at( i ) ); }
        if( i / 3 == 2 ) { side_atom_coord.push_back( coords_xyz.at( i ) ); }
    }

    std::vector< std::vector<double> > reference_frame =
        create_ref_frame( mid_atom_coord, up_atom_coord, side_atom_coord );

    for( int i = 0; i < reference_frame.size(); i++ ) {
        printf( "%.3f %.3f %.3f\n",
                reference_frame[i][0],
                reference_frame[i][1],
                reference_frame[i][2] );
    }

    return 0;
}
