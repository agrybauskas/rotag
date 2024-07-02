%option noyywrap

%{
    #include <cstdlib>
    #include "SelectionParser.h"

    extern int yyparse();
%}

digits [0-9]+

%%

{digits} {
             yyval.integer = atoi(yytext);
             return DIGIT;
         }

%%
