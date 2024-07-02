%define api.prefix select_

%{
    #include <cmath>
    #include <iostream>

    extern int select_lex();

    extern void select_error(char const* msg);
%}

%union
{
    int integer;
}

%token<integer> DIGIT

%start expr

%%

expr: /* empty */
    | DIGIT { std::cout << $1 << std::endl; }
    ;

%%

void select_error(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}
