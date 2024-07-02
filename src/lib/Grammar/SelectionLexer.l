%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>
    #include "SelectionParser.h"

    extern int select_parse();
%}

digits [0-9]+

%%

{digits} {
             select_lval.integer = atoi(select_text);
             return DIGIT;
         }

%%
