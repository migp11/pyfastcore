# -*- coding: utf-8 -*-

from cobra.flux_analysis import single_reaction_deletion
from cobra.core import Reaction, Metabolite, Model
from cobra.test import create_test_model
from fastcore.fastcore import Fastcore


def create_tiny_net():
    model = Model("Tiny_Net")
    metabolites = ["A", "B", "C", "D"]
    net_dict = {"v1": {"A": 2},
                "v2": {"A": -1, "B": 1},
                "v3": {"A": -1, "D": 1},
                "v4": {"A": -1, "C": 1},
                "v5": {"C": -1, "D": 1},
                "v6": {"D": -1}
                }

    bounds = {"v1": [0, 10],
              "v2": [-10, 10],
              "v3": [0, 10],
              "v4": [0, 10],
              "v5": [0, 10],
              "v6": [0, 10]
              }

    to_add = []
    for m in metabolites:
        met = Metabolite(m)
        to_add.append(met)

    model.add_metabolites(to_add)
    to_add = []
    for r in net_dict.keys():
        rxn = Reaction(r)
        rxn.lower_bound = bounds[r][0]
        rxn.upper_bound = bounds[r][1]
        rxn.add_metabolites({model.metabolites.get_by_id(m): s
                             for m, s in net_dict[r].items()})
        to_add.append(rxn)

    model.add_reactions(to_add)

    return model

def test_tiny_net_model():
    model = create_tiny_net()
    for r in model.reactions:
        print(r.id, r.reaction)

    core_reactions = ["v6"]
    fc_builder = Fastcore(model, core_reactions,
                          debug_mode=True)

    fc_builder.fast_core()
    core = fc_builder.get_context_specific_set()
    print("Consistent reaction set size", len(core))
    print("Context specific core:")
    cs_set = fc_builder.get_context_specific_set()
    print(len(cs_set))
    print(cs_set)

def test_textbook_model():
    model = create_test_model('textbook')
    model.reactions.ATPM.lower_bound = 0
    core_reactions = ['Biomass_Ecoli_core']
    penalties = {r: 0 for r in core_reactions}
    for r in model.exchanges:
        penalties[r.id] = 0

    fc_builder = Fastcore(model, core_reactions,
                          penalties=penalties,
                          default_penalty=10,
                          debug_mode=True)

    fc_builder.fast_core()
    core = fc_builder.get_context_specific_set()
    print("Consistent reaction set size", len(core))
    print("Context specific core:")
    cs_set = fc_builder.get_context_specific_set()
    cs_rxn_set = [r for r in cs_set if not r.startswith('EX')]
    print(len(cs_rxn_set) - 1)
    cs_model = fc_builder.build_context_specific_model()
    result = single_reaction_deletion(cs_model)
    print(result)

def test_fastcc():
    model = create_test_model('textbook')

    model.reactions.ATPM.lower_bound = 0
    model.reactions.EX_glc__D_e.lower_bound = -10000
    model.reactions.EX_o2_e.lower_bound = -10000
    blocked = Fastcore.fast_find_blocked(model)
    print(blocked)
    from cobra.flux_analysis import find_blocked_reactions
    blocked = find_blocked_reactions(model)
    print(blocked)

def main():
    test_textbook_model()
    test_fastcc()
    test_tiny_net_model()

main()
