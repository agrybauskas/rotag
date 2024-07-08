%define api.prefix {select_}
%parse-param {Parameters& parameters} {AtomSite& atom_site} {std::set<int64_t>& atom_ids} {int64_t seed} {int64_t group_id}

%code requires{
    #include <string>
    #include <set>

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

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

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

    void set_lex_input(const char*);
    void end_lex_scan();

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
}

%token<num> NUM
%token<dbl> DOUBLE
%token<str> STR SEP ALL MAINCHAIN SIDECHAIN HETATOMS

%start cmd

%%

cmd: /* empty */
    | cmd SEP expr
    | expr
    ;

expr:
    | NUM       { std::printf("%li\n", $1); }
    | DOUBLE    { std::printf("%f\n", $1); }
    | STR       { std::printf("%s\n", $1); }
    | ALL       {
                    std::map<int64_t, Atom> atoms = atom_site.atoms();
                    std::map<int64_t, Atom>::iterator atom_it;
                    for (atom_it = atoms.begin();
                         atom_it != atoms.end();
                         ++atom_it) {
                        atom_ids.insert(atom_it->first);
                    }
                }
    | MAINCHAIN {
                }
    | SIDECHAIN {
                }
    | HETATOMS  {
                    Selector selector;
                    selector["_atom_site.group_pdb"]["HETATM"] = true;
                    std::vector<PDBXVALUE> hetatom_ids =
                        filter(atom_site, selector).ids();
                    for (PDBXVALUE& hetatom_id : hetatom_ids) {
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
