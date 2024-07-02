%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>

    #include "SelectionParser.h"
%}

digits [0-9]+

%%

{digits} {
             select_lval.integer = atoi(select_text);
             return DIGIT;
         }

%%
