%require "3.2"

%define api.prefix {select_}

%parse-param {Parameters& parameters} {AtomSite& atom_site} {std::set<int64_t>& atom_ids} {int64_t seed} {int64_t group_id}

%code requires{
    #include <string>
    #include <set>

    #include "../AtomSite.h"
    #include "../ForceField/Parameters.h"

    struct AtomIDs {
        bool negation = false;
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
    std::vector<int64_t> num_list;
    AtomIDs* atom_ids;
}

%start cmd

%type<atom_ids> cmd expr
%token<atom_ids> MODEL RANGE COMMA NOT LEFT_P RIGHT_P AND OR ALL MAINCHAIN SIDECHAIN HETATOMS
%token<str> NUM DOUBLE STR SEP
%token<num_list> num_oper;

%%

cmd:
    | cmd SEP cmd
    | expr
        {
            for (int64_t atom_id : $1->list) {
                atom_ids.emplace(atom_id);
            }
            delete $1;
        }
    ;

expr:
    | ALL
        {
            $$ = new AtomIDs();
            for (PDBXVALUE& atom_id : atom_site.ids()) {
                $$->list.emplace((int64_t) atom_id);
            }
        }
    | MAINCHAIN
        {
            $$ = new AtomIDs();
            std::vector<PDBXVALUE> mainchain_atom_ids =
                filter(atom_site, parameters.MAINCHAIN_ATOM_NAMES).ids();
            for (PDBXVALUE& mainchain_atom_id : mainchain_atom_ids) {
                $$->list.emplace((int64_t) mainchain_atom_id);
            }
        }
    | SIDECHAIN
        {
            $$ = new AtomIDs();
            std::vector<PDBXVALUE> sidechain_atom_ids =
                filter(atom_site, parameters.SIDECHAIN_ATOM_NAMES).ids();
            for (PDBXVALUE& sidechain_atom_id : sidechain_atom_ids) {
                $$->list.emplace((int64_t) sidechain_atom_id);
            }
        }
    | HETATOMS
        {
            $$ = new AtomIDs();
            Selector selector = {{{"_atom_site.group_pdb", {"HETATM"}}}};
            std::vector<PDBXVALUE> hetatom_ids =
                filter(atom_site, selector).ids();
            for (PDBXVALUE& hetatom_id : hetatom_ids) {
                $$->list.emplace((int64_t) hetatom_id);
            }
        }
    | expr AND expr
        {
            $$ = new AtomIDs();
            std::set<int64_t> atom_ids_1 = $1->list;
            for (int64_t atom_id_1 : atom_ids_1) {
                std::set<int64_t> atom_ids_2 = $3->list;
                if (auto search = atom_ids_2.find(atom_id_1);
                    search != atom_ids_2.end()) {
                    $$->list.emplace(atom_id_1);
                }
            }
        }
    | expr OR expr
        {
            $$ = new AtomIDs();
            std::set<int64_t> atom_ids_1 = $1->list;
            for (int64_t atom_id_1 : atom_ids_1) {
                $$->list.emplace(atom_id_1);
            }
            std::set<int64_t> atom_ids_2 = $3->list;
            for (int64_t atom_id_2 : atom_ids_2) {
                $$->list.emplace(atom_id_2);
            }
        }
    | LEFT_P expr RIGHT_P
        {
            $$ = new AtomIDs();
            for (int64_t atom_id : $2->list) {
                $$->list.emplace(atom_id);
            }
        }
    | NOT expr
        {
            $$ = new AtomIDs();
            $2->negation = true;
            for (int64_t atom_id : $2->list) {
                $$->list.emplace(atom_id);
            }
        }
    | MODEL num_oper
        {
            $$ = new AtomIDs();
            Selector selector;
            //for () {
            //
            //}
            //for (PDBXVALUE& atom_id : atom_ids) {
                //$$->list.emplace((int64_t) atom_id);
            //}
        }
    ;

//any_oper:
//    | any_oper COMMA any_oper
//    | num_oper
//    | str_oper
//    ;

num_oper:
    | num_oper COMMA num_oper
    | num_oper RANGE num_oper
    | NUM
    ;

//str_oper:
//    | str_oper COMMA str_oper
//    | STR
//    ;

    /* | DOUBLE    {} */
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
