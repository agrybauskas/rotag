%{
#include "SelectionParser.h"
%}
%%
[0-9]+ {return NUM;}
[A-Za-z]+ {return STR;}
%%