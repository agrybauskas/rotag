%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>

    #include "SelectionParser.h"
%}

num   [0-9]+
space [ \t\n]+
separator [\;]

%%

{num}       {
                select_lval.num = atoi(select_text);
                return NUM;
            }
{space}     { return SPACE; }
{separator} { return SEPARATOR; }

%%

void set_lex_input(const char* input) {
    yy_scan_string(input);
}

void end_lex_scan(void) {
    yy_delete_buffer(YY_CURRENT_BUFFER);
}
