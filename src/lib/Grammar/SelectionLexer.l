%{
#include "SelectionParser.h"
%}
%option noyywrap
%%
[0-9]+ {return NUM;}
[A-Za-z]+ {return STR;}
%%
