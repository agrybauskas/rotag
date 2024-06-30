%defines

%{
    #include <string>

    extern int yylex();

    extern void yyerror(char const* msg);
%}

/* %token<std::string> COMMA */
%token<int64_t> NUM
/* %token<std::string> STR */
/* %left COMMA */
/* %left NUM STR */
/* %left ".." "=" */
/* %left 'around' 'rand' 'angles' */
/* %right 'byres' 'expand' */
/* %left '!' */
/* %left '||' '&&' */
/* %left ':' */

%%

exp:
exp:
    | NUM {std::cout << $1 << std::endl;}
    ;

/* exp */
/*     : any_ope */
/* ; */

/* any_ope */
/*     /\* : any_ope COMMA any_ope *\/ */
/*     : num_ope */
/*     | str_ope */
/* ; */

/* num_ope */
/*     : NUM */
/* ; */

/* str_ope */
/*     : STR */
/* ; */

%%

void yyerror(char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

/* void selection_parser(AtomSite& atom_site, std::string query) { */
/* } */
