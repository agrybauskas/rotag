%option noyywrap
%option prefix="select_"

%{
    #include "SelectionParser.h"
%}

all       (all|\*)
mainchain (mainchain|mc)
sidechain (sidechain|sc)
hetatoms  (hetatoms|het)
num       [0-9]+
double    [-+]?([0-9]*[\.])?[0-9]+([eE][-+]?[0-9]+)?
str       [A-Za-z\.]+
space     [ \t]+
separator [\;]

%%

{all}       { return ALL; }
{mainchain} { return MAINCHAIN; }
{sidechain} { return SIDECHAIN; }
{hetatoms}  { return HETATOMS; }
{num}       {
                ;; select_lval.num = std::stoll(select_text);
                return NUM;
            }
{double}    {
                ;; select_lval.dbl = std::stod(select_text);
                return DOUBLE;
            }
{str}       {
                ;; select_lval.str = select_text;
                return STR;
            }
{space}     { /* ignore */ }
{separator} { return SEP; }

%%

void set_lex_input(const char* input) {
    yy_scan_string(input);
}

void end_lex_scan() {
    yy_delete_buffer(YY_CURRENT_BUFFER);
    yylex_destroy();
}
