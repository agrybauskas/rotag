%define api.prefix select_

%{
    #include <cmath>
    #include <iostream>
    #include <string>
    #include <vector>

    extern int select_lex();

    extern void select_error(char const* msg);
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

std::vector<int64_t> select_parse(std::string cmd) {
    return std::vector<int64_t>({1, 2, 3});
}

void select_error(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}
