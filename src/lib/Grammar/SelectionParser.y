%code requires{
  #include <string>
  #include <vector>

  #include "../PDBxParser.h"

  void selection_parser(AtomSite atom_site, std::string query);
}

%{
#include <stdio.h>
#include <stdlib.h>

#include "SelectionLexer.h"

int yyparse();
int yylex(void);
void yyerror(std::string message){}
%}

%token NUM STR COMMA

/* %token NUM FLOAT STR */
%left COMMA
%left NUM STR
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
    : any_ope COMMA any_ope
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
