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
%token<str> STR SEP

%start cmd

%%

cmd: /* empty */
    | cmd SEP expr
    | expr
    ;

expr:
    | NUM    { std::printf( "%li\n", $1 ); }
    | DOUBLE { std::printf( "%f\n", $1 ); }
    | STR    { std::printf( "%s\n", $1 ); }
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
