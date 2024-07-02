%option noyywrap

%{
    #define YY_DECL int select_lex()

    #include <cstdlib>
    #include "SelectionParser.h"

    extern int yyparse();
%}

digits [0-9]+

%%

{digits} {
             yylval.integer = atoi(yytext);
             return DIGIT;
         }

%%
