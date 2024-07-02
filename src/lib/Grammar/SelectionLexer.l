%option noyywrap
%option prefix="select_"

%{
    #include "SelectionParser.h"
%}

all       (all|\*)
num       [0-9]+
double    [-+]?([0-9]*[\.])?[0-9]+([eE][-+]?[0-9]+)?
str       [A-Za-z\.]+
space     [ \t]+
separator [\;]

%%

{all}       { return ALL; }
{num}       {
                select_lval.num = std::stoll(select_text);
                return NUM;
            }
{double}    {
                select_lval.dbl = std::stod(select_text);
                return DOUBLE;
            }
{str}       {
                select_lval.str = select_text;
                return STR;
            }
{space}     { /* ignore */ }
{separator} { return SEP; }

%%

void set_lex_input(const char* input) {
    yy_scan_string(input);
}

void end_lex_scan(void) {
    yy_delete_buffer(YY_CURRENT_BUFFER);
}
