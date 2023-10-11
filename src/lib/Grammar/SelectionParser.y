%{
#include "SelectionLexer.h"

#include <stdio.h>
#include <stdlib.h>

void yyerror(const char *msg);
%}

%token NUM
%token FLOAT
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
    : any_ope "," any_ope
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
void yyerror(const char *msg) {}
