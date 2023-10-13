%require "3.2"

%{
#include "SelectionLexer.h"

#include <stdio.h>
#include <stdlib.h>
%}

%code requires{
#include <string>
#include <vector>

#include "../PDBxParser.h"

std::vector<std::string>
selection_parser(AtomSite atom_site, std::string cmd_line);
}

%token NUM
%token FLOAT
%token STR

%left ","
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

std::vector<std::string>
selection_parser(AtomSite atom_site, std::string cmd_line) {
  return std::vector<std::string>{};
}
