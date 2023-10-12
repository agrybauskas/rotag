%{
#include "SelectionParser.h"
%}

%option noyywrap yylineno

%%

[0-9]+    {return NUM;}
[A-Za-z]+ {return STR;}

%%
