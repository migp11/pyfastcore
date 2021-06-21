import json

def read_json(json_fname):
    result = None
    with open(json_fname) as fh:
        result = json.load(fh)
    return result


def set_medium(model, medium, skip_absent=True, inplace=False):
    if not inplace:
        model = model.copy()

    def set_active_bound(reaction, bound):
        if reaction.reactants:
            reaction.lower_bound = -bound
        elif reaction.products:
            reaction.upper_bound = bound

    # Set the given media bounds
    media_rxns = list()
    for rxn_id, bound in medium.items():
        if rxn_id not in model.reactions and skip_absent:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        media_rxns.append(rxn)
        set_active_bound(rxn, bound)

    boundary_rxns = set([r for r in model.exchanges if r.id.startswith('EX')])
    media_rxns = set(media_rxns)

    # Turn off reactions not present in media
    for rxn in (boundary_rxns - media_rxns):
        set_active_bound(rxn, 0)

    if not inplace:
        return model
