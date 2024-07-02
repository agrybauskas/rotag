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

%token<integer> DIGIT

%start expr

%%

expr: /* empty */
    | DIGIT { std::cout << $1 << std::endl; }
    ;

%%

void yyerror(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}
