%define api.prefix select_
%parse-param {std::vector<int64_t>& atom_ids}

%code requires{
    #include <vector>

    #include "../PDBxParser.h"

    std::vector<int64_t> selection_parser(AtomSite&, std::string);
}

%{
    #include <cmath>
    #include <iostream>
    #include <string>
    #include <vector>

    extern int select_lex();

    extern void select_error(std::vector<int64_t>&, char const*);
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

void select_error(std::vector<int64_t>& atom_ids, char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

std::vector<int64_t> selection_parser(AtomSite& atom_site, std::string cmd) {
    std::vector<int64_t> atom_ids = {1, 2, 3};
    return atom_ids;
}
