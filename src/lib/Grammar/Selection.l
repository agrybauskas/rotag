;; %option yylineno
%option noyywrap

%{
    #include "SelectionParser.h"

    extern int yyparse();
%}

;; COMMA ","
NUM [0-9]+
;; STR [A-Za-z]+

%%

{COMMA}

%%
