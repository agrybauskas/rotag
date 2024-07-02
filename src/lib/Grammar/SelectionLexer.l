%option noyywrap
%option prefix="select_"

%{
    #include <cstdlib>

    #include "SelectionParser.h"
%}

NUM [0-9]+

%%

{NUM} {
          select_lval.num = atoi(select_text);
          return NUM;
      }

%%

void set_lex_input(const char* input) {
    yy_scan_string(input);
}

void end_lex_scan(void) {
    yy_delete_buffer(YY_CURRENT_BUFFER);
}
