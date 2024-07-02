%{
    #include <cmath>
    #include <iostream>

    extern int yylex();

    extern void yyerror(char const* msg);
%}

%union
{
    int integer;
}

%token<int> DIGIT

%start expr

%%

exp: /* empty */
    | DIGIT {std::cout << $1 << std::endl;}
    ;

%%

void yyerror(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}
