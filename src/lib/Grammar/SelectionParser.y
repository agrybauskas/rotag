%define api.prefix {select_}
%parse-param {Parameters& parameters} {AtomSite& atom_site} {std::set<int64_t>& atom_ids} {int64_t seed} {int64_t group_id}

%code requires{
    #include <string>
    #include <set>
    #include <vector>

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

    class Set;

    std::set<int64_t> selection_parser(Parameters&,
                                       AtomSite&,
                                       std::string,
                                       int64_t,
                                       int64_t);
}

%{
    #include <iostream>
    #include <string>
    #include <set>
    #include <vector>

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

    void set_lex_input(const char*);
    void end_lex_scan();

    class Set {
     public:
        std::set<int64_t> data;
        Set() { data = {}; }
        void add(int64_t atom_id) {
            this->data.insert(atom_id);
        }
    };

    extern int select_lex();

    extern void select_error(Parameters&,
                             AtomSite&,
                             std::set<int64_t>&,
                             int64_t,
                             int64_t,
                             char const*);
%}

%union
{
    int64_t num;
    char* str;
    double dbl;
    Set* set;
}

%start cmd

%token<num> NUM
%token<dbl> DOUBLE
%token<str> STR SEP ALL MAINCHAIN SIDECHAIN HETATOMS

%%

cmd:
    | cmd SEP expr
    | expr
    ;

expr:
    | NUM       { std::printf("%li\n", $1); }
    | DOUBLE    { std::printf("%f\n", $1); }
    | STR       { std::printf("%s\n", $1); }
    | ALL       {
                    for (PDBXVALUE& id : atom_site.ids()) {
                        /* $$->add((int64_t) id); */
                        atom_ids.insert((int64_t) id);
                    }
                }
    | MAINCHAIN {
                }
    | SIDECHAIN {
                }
    | HETATOMS  {
                    Selector selector =
                        {{"_atom_site.group_pdb", {{"HETATM", true}}}};
                    std::vector<PDBXVALUE> hetatom_ids =
                        filter(atom_site, selector).ids();
                    for (PDBXVALUE& hetatom_id : hetatom_ids) {
                        /* $$->add((int64_t) hetatom_id); */
                        atom_ids.insert((int64_t) hetatom_id);
                    }
                }
    ;

%%

void select_error(Parameters& parameters,
                  AtomSite&,
                  std::set<int64_t>& atom_ids,
                  int64_t seed,
                  int64_t group_id,
                  char const* msg) {
    std::cout << "Syntax Error: " << msg << std::endl;
}

std::set<int64_t> selection_parser(Parameters& parameters,
                                   AtomSite& atom_site,
                                   std::string cmd,
                                   int64_t seed,
                                   int64_t group_id) {
    std::set<int64_t> atom_ids = {};
    set_lex_input(cmd.c_str());
    int retval = select_parse(parameters, atom_site, atom_ids, seed, group_id);
    end_lex_scan();
    return atom_ids;
}
