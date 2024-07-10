%require "3.2"

%define api.prefix {select_}

%parse-param {Parameters& parameters} {AtomSite& atom_site} {std::set<int64_t>& atom_ids} {int64_t seed} {int64_t group_id}

%code requires{
    #include <string>
    #include <set>

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

    struct Data {
        std::set<int64_t> list;
    };

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
    char* str;
    Data* data;
}

%start cmd

%token<data> AND ALL MAINCHAIN SIDECHAIN HETATOMS
%token<str> NUM DOUBLE STR SEP
%type<data> cmd expr

%%

cmd:
    | cmd SEP cmd
    | expr
        {
            $$ = $1;
            for (int64_t atom_id : $$->list) {
                atom_ids.emplace(atom_id);
            }
            delete $$;
        }
    ;

expr:
    | ALL
        {
            for (PDBXVALUE& atom_id : atom_site.ids()) {
                $$->list.emplace((int64_t) atom_id);
            }
        }
    | MAINCHAIN
        {

        }
    | SIDECHAIN
        {

        }
    | HETATOMS
        {
            Selector selector =
                {{"_atom_site.group_pdb", {{"HETATM", true}}}};
            std::vector<PDBXVALUE> hetatom_ids =
                filter(atom_site, selector).ids();
            for (PDBXVALUE& hetatom_id : hetatom_ids) {
                $$->list.emplace((int64_t) hetatom_id);
            }
        }
    | expr AND expr
        {
            $$ = $2;
            std::set<int64_t> atom_ids_1 = $1->list;
            for (int64_t atom_id_1 : atom_ids_1) {
                std::set<int64_t> atom_ids_2 = $3->list;
                if (auto search = atom_ids_2.find(atom_id_1);
                    search != atom_ids_2.end()) {
                    $$->list.emplace(atom_id_1);
                }
            }
        }
    /* | NUM       {} */
    /* | DOUBLE    {} */
    /* | STR       {} */
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
