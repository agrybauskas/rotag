%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>

    #include "SelectionParser.h"
%}

NUM [0-9]+

%%

{NUM} {
          select_lval.num = atoi(select_text);
          return NUM;
      }

%%
