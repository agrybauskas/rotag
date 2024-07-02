%define api.prefix select_
%parse-param {AtomSite& atom_site} {std::vector<int64_t>& atom_ids}

%code requires{
    #include <string>
    #include <vector>

    #include "../PDBxParser.h"

    std::vector<int64_t> selection_parser(AtomSite&, std::string);
}

%{
    #include <iostream>
    #include <string>
    #include <vector>

    #include "../PDBxParser.h"

    void set_lex_input(const char* input);
    void end_lex_scan(void);

    extern int select_lex();

    extern void select_error(AtomSite&, std::vector<int64_t>&, char const*);
%}

%union
{
    int64_t num;
    char* str;
    double dbl;
}

%token<num> NUM
%token<dbl> DOUBLE
%token<str> STR SPACE SEP

%start cmd

%%

cmd: /* empty */
    | cmd SEP expr
    | expr
    ;

expr:
    /* | NUM SPACE { std::cout << $1 << std::endl; } */
    /* | SPACE NUM { std::cout << $2 << std::endl; } */
    /* | SPACE NUM SPACE { std::cout << $2 << std::endl; } */
    /* | STR SPACE { std::cout << $1 << std::endl; } */
    /* | SPACE STR { std::cout << $2 << std::endl; } */
    /* | SPACE STR SPACE { std::cout << $2 << std::endl; } */
    | NUM    { std::cout << $1 << std::endl; }
    | DOUBLE { std::cout << $1 << std::endl; }
    | STR    { std::cout << $1 << std::endl; }
    ;

%%

void select_error(AtomSite&, std::vector<int64_t>& atom_ids, char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

std::vector<int64_t> selection_parser(AtomSite& atom_site, std::string cmd) {
    std::vector<int64_t> atom_ids = {};
    set_lex_input(cmd.c_str());
    int retval = select_parse(atom_site, atom_ids);
    end_lex_scan();
    return atom_ids;
}
