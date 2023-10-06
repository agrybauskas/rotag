/* SelectionParser.l */

%{
#include "SelectionParser.h"
%}

%token NUM
%left ","
/* %left '..' '=' */
/* %left 'around' 'rand' 'angles' */
/* %right 'byres' 'expand' */
/* %left '!' */
/* %left '||' '&&' */
/* %left ':' */

%%
line
    : NUM
%%
