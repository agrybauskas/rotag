%define api.prefix select_
%parse-param {std::vector<int>& atom_ids}

%{
    #include <cmath>
    #include <iostream>
    #include <string>
    #include <vector>

    #include "../PDBxParser.h"

    extern int select_lex();

    extern void select_error(std::vector<int>&, char const*);

    std::vector<int> selection_parser(AtomSite&, std::string);
%}

%union
{
    int integer;
}

%token<integer> DIGIT

%start cmd

%%

cmd: /* empty */
    | DIGIT { std::cout << $1 << std::endl; }
    ;

%%

void select_error(std::vector<int>& atom_ids, char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

std::vector<int> selection_parser(AtomSite& atom_site, std::string cmd) {

}
