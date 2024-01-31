%{
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <vector>

#include "SelectionLexer.h"

#include "../PDBxParser.h"

void selection_parser(AtomSite atom_site, std::string query);

extern "C" {
  int yylex(void);
  int yyerror(char *message);
}
%}

%token NUM STR

/* %token NUM FLOAT STR */
%left COMMA
/* %left ".." "=" */
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

void selection_parser(AtomSite atom_site, std::string query) {
}
