%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>
    #include <string>
    #include <vector>

    #include "SelectionParser.h"

    extern std::vector<int64_t> select_parse(std::string);
%}

digits [0-9]+

%%

{digits} {
             select_lval.integer = atoi(select_text);
             return DIGIT;
         }

%%
