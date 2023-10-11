/* SelectionParser.l */

%{
#include "SelectionParser.h"
%}

%token NUM
%token STR

%left ","
%left ".." "="
/* %left 'around' 'rand' 'angles' */
/* %right 'byres' 'expand' */
/* %left '!' */
/* %left '||' '&&' */
/* %left ':' */

%%
exp
    : any_ope
;

any_ope
    : num_ope "," str_ope
    | str_ope "," num_ope
    | num_ope
    | str_ope
;

num_ope
    : NUM
;

str_ope
    : STR
;

%%
