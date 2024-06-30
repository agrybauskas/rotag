%defines

%{
    #include <cmath>
    #include <iostream>
    #include <string>

    extern int yylex();

    extern void yyerror(char const* msg);
%}

%token<int> NUM

%%

exp: /* empty */
    | NUM {std::cout << $1 << std::endl;}
    ;

%%

void yyerror(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

/* void selection_parser(AtomSite& atom_site, std::string query) { */
/* } */
