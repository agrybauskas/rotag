%option noyywrap

%{
    #include <cstdlib>
    #include "SelectionParser.h"

    extern int yyparse();
%}

digits [0-9]+

%%

{digits} {
             yyval = atoi(yytext);
             return NUM;
         }

%%
