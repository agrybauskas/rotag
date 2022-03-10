#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/regex.hpp>
#include "PDBxParser.h"

/*
  Obtains pdbx loops for a specified categories.
  Input:
      pdbx_file - PDBx file path;
      categories - list of specified categories.
      read_until_end - boolean flag for reading whole pdbx file or stdin.
  Output:
      pdbx_loop_data_list - data structure for loop data or list of data structure.
*/

std::vector<pdbx_loop_data> obtain_pdbx_loop( std::string pdbx_file,
                                              std::vector<std::string> categories,
                                              bool read_until_end = false )
{
    std::vector< std::vector<std::string> > current_categories;
    std::vector< std::vector< std::vector<std::string> > > attributes;
    /* Will be used for storing atom data temporarily.*/
    std::vector< std::vector< std::vector<std::string> > > data;

    boost::regex category_regex(
        "(" + boost::algorithm::join( categories , "|") + ")[.](.+)\n?$"
    );
    boost::regex data_regex( "^data_" );
    boost::regex reading_end_regex( "^_|loop_|#" );

    bool is_reading_lines = false; /* Starts/stops reading lines at certain flags. */

    std::ifstream fh( pdbx_file );
    std::string line;
    boost::smatch match;
    if( fh.is_open() ) {
        while( getline( fh, line ) ) {
            if( boost::regex_search( line, match, data_regex ) ||
                current_categories.empty() ) {
                current_categories.push_back( std::vector<std::string>() );
                attributes.push_back( std::vector< std::vector<std::string> >() );
                data.push_back( std::vector< std::vector<std::string> >() );
            } else if( boost::regex_search( line, match, category_regex ) ) {
                if( current_categories.back().empty() ||
                    current_categories.back().back() != match[1] ) {
                    std::string category( match[1] );
                    current_categories.back().push_back( category );
                    attributes.back().push_back( std::vector<std::string>() );
                    data.back().push_back( std::vector<std::string>() );
                }
                std::string attribute( match[2] );
                boost::trim( attribute );
                attributes.back().back().push_back( attribute );
                is_reading_lines = true;
            } else if( is_reading_lines &&
                       boost::regex_search( line, match, reading_end_regex ) ) {
                if( current_categories.size() == current_categories.size() &&
                    ! read_until_end ) {
                    break;
                }
                is_reading_lines = false;
            } else if( is_reading_lines ) {
                std::vector<std::string> data_split;
                boost::split( data_split, line, boost::is_any_of(" ")  );
                for( int i = 0; i < data_split.size(); i++ ) {
                    if( data_split[i] != "" ) {
                        data.back().back().push_back( data_split[i] );
                    }
                }
            }
        }
        fh.close();
    }

    /* Generates hash from three lists. */
    std::vector<pdbx_loop_data> pdbx_loop_data_list;
    for( int i = 0; i < current_categories.size(); i++ ) {
        pdbx_loop_data pdbx_loop_data;
        for( int j = 0; j < current_categories.back().size(); j++ ) {
            pdbx_loop_data[current_categories[i][j]]["attributes"] = attributes[i][j];
            pdbx_loop_data[current_categories[i][j]]["data"] = data[i][j];
        }

        pdbx_loop_data_list.push_back( pdbx_loop_data );
    }

    return pdbx_loop_data_list;
}

/*
  Takes PDBx loop and converts it to hash of hashes where the first key is
  unique.
  Input:
      pdbx_loop_data - data structure (from obtain_pdbx_loop);
      unique_keys - combination of attribute data that serves as unique key.
      read_until_end - reads whole pdbx file or stdin.
      Ex.: [ 'id' ]
  Output:
      pdbx_loop_unique - special data structure.
*/

pdbx_unique_entry pdbx_loop_unique( std::vector<pdbx_loop_data> pdbx_loop_data,
                                    std::vector<std::string> unique_keys,
                                    bool read_until_end = false )
{

}
