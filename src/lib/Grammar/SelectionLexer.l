%{
#include "SelectionParser.h"
%}

%option yylineno
%option noyywrap

%%

","       {return COMMA;}
[0-9]+    {return NUM;}
[A-Za-z]+ {return;}

%%
