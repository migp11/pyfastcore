# -*- coding: utf-8 -*-

"""This is a Python implementation of the fastcore algorithm [1]
[1] Vlassis, N., Pacheco, M. P., & Sauter, T. (2014). PLoS Computational Biology, 10(1)
<http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424>

"""
import six
from sympy import S


class Fastcore(object):

    AUXILIARY_VAR_PREFIX = "Z_var_"
    AUXILIARY_CONS_PREFIX = "Z_const_"

    def __init__(self, model, core_reactions, penalties=None, epsilon=1e-4, scaling_factor=1e5,
                 default_penalty=10.0, check_consistency=True, debug_mode=False, tolerance=1e-7):

        """
        :param model: A cobrapy model used as the universal network for the reconstruction.
        :param core_reactions: A list of reaction that should be included in the context-specific model
        :param penalties: A dictionary including the cost of including a non-core reactions
        :param epsilon: A scalar value used as a lower bound to force the flux through the core_reactions
        :param scaling_factor: A scalar value used to re-scale epsilon in LP10 (see original publication)
        :param default_penalty: The default cost value of including non-core reaction.
        :param consistency_check: If True, check model consistency and prune blocked reaction.
        :param debug_mode: If True, print debugging information.
        :param tolerance: A scalar value used as the zero cutoff
        """

        self._original_model = model
        self._model = self._original_model.copy()

        self.all_reactions = frozenset({r.id for r in model.reactions})

        assert len(core_reactions) > 0

        if hasattr(core_reactions[0], 'id'):
            core_reactions = [r.id for r in core_reactions]

        for r in core_reactions:
            if r in self.all_reactions:
                continue
            raise Exception("Reaction %s not included in model %s" % (r, model.id))

        self.core_reactions = frozenset(core_reactions)

        self._default_penalty = default_penalty
        self.penalties = {}
        self._initialize_penalties(penalties)

        self._tolerance = tolerance
        self._epsilon = epsilon
        self._scaling_factor = scaling_factor
        self.LP7 = None
        self.LP10 = None
        self._debug_mode = debug_mode
        self.bounds = dict()
        for r in self._model.reactions:
            self.bounds[r.id] = r.bounds

        self.blocked_reactions = frozenset()
        self.gap_metabolites = frozenset()
        self._consistent_sets = None

        print("===========================================================")
        print("Initializing Fastcore Builder using")
        print("Model: %s" % model.id)
        print("- Nº of reactions: %i" % len(self.all_reactions))
        print("- Nº of core reactions: %i" % len(self.core_reactions))

        if check_consistency:
            print("Checking network consistency (may take some minutes)")
            self._consistency_check()
            if not self.is_consistent:
                self._prune_inconsistencies()
                print("Warning, original model contains inconsistent reactions")
                print("- Nº of blocked reactions: %i" % len(self.blocked_reactions))
                print("- Nº of gap metabolites: %i" % len(self.gap_metabolites))
                print("\nPruning model inconsistencies")
                print("- Nº of reactions after pruning: %i" % len(self.consistent_reactions))
                print("- Nº of core reactions after pruning: %i\n" % len(self.consistent_core_reactions))

        print("Initializing LP7 and LP10")
        self.LP7 = Fastcore.create_optlang_lp7(self._model,
                                               the_reactions=self.consistent_core_reactions,
                                               epsilon=epsilon
                                               )

        self.LP10 = Fastcore.create_optlang_lp10(self._model,
                                                 self.penalties,
                                                 scaling_factor=scaling_factor)
        print("Fastcore builder ready!")
        print("===========================================================")

    @property
    def tolerance(self):
        return self._tolerance

    @property
    def epsilon(self):
        return self._epsilon

    @property
    def scaling_factor(self):
        return self._scaling_factor

    @property
    def blocked_core_reactions(self):
        return self.core_reactions & self.blocked_reactions

    @property
    def consistent_reactions(self):
        return self.all_reactions - self.blocked_reactions

    @property
    def consistent_core_reactions(self):
        return self.consistent_reactions & self.core_reactions

    @property
    def is_consistent(self):
        return self.blocked_core_reactions == 0

    @property
    def consistent_non_core_reactions(self):
        return self.consistent_reactions - self.core_reactions

    # Class methods

    @staticmethod
    def create_reactions_dict(cobra_model):
        rxns_dict = {}
        for r in cobra_model.reactions:
            rxns_dict[r.id] = {m.id: c for m, c in r.metabolites.items()}
        return rxns_dict

    @staticmethod
    def flip_variable(rxn_id, lp_model, reaction_dict):
        var = lp_model.variables.get(rxn_id)
        lb, ub = var.lb, var.ub
        var.ub = -lb
        var.lb = -ub
        for cont_id, coeff in reaction_dict.items():
            const = lp_model.constraints.get(cont_id)
            const.set_linear_coefficients({var: coeff})
        lp_model.update()

    @staticmethod
    def create_optlang_fba(cobra_model):
        optlang_interface = cobra_model.solver.interface
        lp_model = optlang_interface.Model()
        lp_model.name = 'FBA'
        to_add = []
        for met in cobra_model.metabolites:
            to_add += [optlang_interface.Constraint(S.Zero, name=met.id, lb=0, ub=0)]

        lp_model.add(to_add)

        constraint_terms = {}
        for reaction in cobra_model.reactions:
            lb, ub = reaction.bounds
            variable = optlang_interface.Variable(reaction.id, lb=lb, ub=ub)
            variable.is_flux = True
            lp_model.add(variable)

            for metabolite in reaction.metabolites:
                constraint = lp_model.constraints[metabolite.id]
                if constraint not in constraint_terms:
                    constraint_terms[constraint] = dict()
                constraint_terms[constraint][variable] = reaction.get_coefficient(metabolite)

        lp_model.update()

        for constraint, terms in six.iteritems(constraint_terms):
            constraint.set_linear_coefficients(terms)

        objective_coefficients = {lp_model.variables.get(r.id): r.objective_coefficient
                                  for r in cobra_model.reactions}

        lp_model.objective.set_linear_coefficients(objective_coefficients)
        lp_model.objective.direction = 'max'
        lp_model.update()

        return lp_model

    @staticmethod
    def create_optlang_lp7(cobra_model, the_reactions=[], epsilon=1e-4):
        optlang_interface = cobra_model.solver.interface
        lp_model = Fastcore.create_optlang_fba(cobra_model)
        lp_model.name = 'Fastcore LP7'

        if len(the_reactions) == 0:
            the_reactions = [r.id for r in cobra_model.reactions]
        else:
            if hasattr(list(the_reactions)[0], 'id'):
                the_reactions = [r.id for r in the_reactions]

        objective_coefficients = {v: 0.0 for v in lp_model.variables}
        to_add = []
        for rxn_id in the_reactions:
            z_var_id = Fastcore.AUXILIARY_VAR_PREFIX + rxn_id
            z_var = optlang_interface.Variable(z_var_id, lb=0, ub=epsilon)
            lp_model.add(z_var)
            objective_coefficients[z_var] = 1.

            v_var = lp_model.variables.get(rxn_id)
            z_constraint_id = Fastcore.AUXILIARY_CONS_PREFIX + rxn_id
            constraint = optlang_interface.Constraint(v_var - z_var, lb=0, name=z_constraint_id)
            to_add.append(constraint)

        lp_model.add(to_add)
        lp_model.objective.set_linear_coefficients(objective_coefficients)
        lp_model.objective.direction = 'max'
        lp_model.update()

        return lp_model

    @staticmethod
    def update_optlang_lp7(lp_model, core_reactions, relaxed_lb=-100000.):

        for c in lp_model.constraints:
            if not c.name.startswith(Fastcore.AUXILIARY_CONS_PREFIX):
                continue
            c.lb = relaxed_lb

        objective_coefficients = {v: 0 for v in lp_model.variables}
        for rxn_id in core_reactions:
            z_var_id = Fastcore.AUXILIARY_VAR_PREFIX + rxn_id
            var = lp_model.variables.get(z_var_id)
            objective_coefficients[var] = 1

            z_constraint_id = Fastcore.AUXILIARY_CONS_PREFIX + rxn_id
            constraint = lp_model.constraints.get(z_constraint_id)
            constraint.lb = 0

        lp_model.objective.set_linear_coefficients(objective_coefficients)
        lp_model.update()

    @staticmethod
    def create_optlang_lp10(cobra_model, penalties, scaling_factor=1e5):
        lp_model = Fastcore.create_optlang_fba(cobra_model)
        lp_model.name = 'Fastcore LP10'

        optlang_interface = cobra_model.solver.interface
        objective_coefficients = {v: 0.0 for v in lp_model.variables}
        to_add = []
        for v_var in lp_model.variables:
            # Not really sure why, but I have to rescale de variables
            # an order of magnitude less than the scaling factor
            v_var.lb *= scaling_factor / 100.
            v_var.ub *= scaling_factor / 100.
            if penalties[v_var.name] == 0:
                continue

            z_var_id = Fastcore.AUXILIARY_VAR_PREFIX + v_var.name
            z_var = optlang_interface.Variable(z_var_id, lb=0, ub=v_var.ub)
            objective_coefficients[z_var] = penalties[v_var.name]
            lp_model.add(z_var)

            constraint_id = "_".join([Fastcore.AUXILIARY_CONS_PREFIX, "LHS", v_var.name])
            lhs_constraint = optlang_interface.Constraint(z_var + v_var, lb=0, name=constraint_id)
            to_add.append(lhs_constraint)

            constraint_id = "_".join([Fastcore.AUXILIARY_CONS_PREFIX, "RHS", v_var.name])
            rhs_constraint = optlang_interface.Constraint(z_var - v_var, lb=0, name=constraint_id)
            to_add.append(rhs_constraint)

        lp_model.add(to_add)
        lp_model.objective.set_linear_coefficients(objective_coefficients)
        lp_model.objective.direction = 'min'
        lp_model.update()

        return lp_model

    @staticmethod
    def update_optlang_lp10(lp_model, core_reactions, penalties, bounds,
                            scaling_factor=1e-5, epsilon=1e-4):

        objective_coefficients = {}
        for var in lp_model.variables:
            if var.name.startswith(Fastcore.AUXILIARY_VAR_PREFIX):
                rxn_var_id = var.name[6:]
                objective_coefficients[var] = penalties[rxn_var_id]
            elif var.name in core_reactions:
                var.lb = epsilon * scaling_factor
            else:
                var.lb = bounds[var.name][0] * scaling_factor / 100.

        lp_model.objective.set_linear_coefficients(objective_coefficients)
        lp_model.update()

    @staticmethod
    def fastcc(cobra_model, reactions_to_check=[], epsilon=1e-4, tolerance=1e-8, debug=False):
        reactions_dict = Fastcore.create_reactions_dict(cobra_model)
        consistent_rxns = set()
        eps = 0.99 * epsilon

        if len(reactions_to_check) == 0:
            reactions_to_check = set([r.id for r in cobra_model.reactions])
        lp_model = Fastcore.create_optlang_lp7(cobra_model, the_reactions=reactions_to_check, epsilon=epsilon)

        irreversible_reactions = {r.id for r in cobra_model.reactions if not r.reversibility}
        to_check = reactions_to_check & irreversible_reactions
        Fastcore.update_optlang_lp7(lp_model, to_check)
        lp_model.optimize()
        if lp_model.status == 'infeasible':
            raise Exception("INFEASIBLE LP7")
        active_reactions = {r for r in reactions_to_check
                            if abs(lp_model.variables.get(r).primal) - eps >= tolerance}

        consistent_rxns |= active_reactions
        reactions_to_check -= consistent_rxns
        reactions_to_check -= irreversible_reactions

        rxns_subset = {}
        flipped = False
        singleton = False
        while reactions_to_check:
            reactions_to_check -= consistent_rxns
            if debug:
                print("Nº of remaining reactions to check:", len(reactions_to_check))

            if singleton:
                rxn = sorted(reactions_to_check)[0]
                Fastcore.update_optlang_lp7(lp_model, {rxn})
            else:
                Fastcore.update_optlang_lp7(lp_model, reactions_to_check)

            lp_model.optimize()
            if lp_model.status == 'infeasible':
                raise Exception("INFEASIBLE LP7")
            active_reactions = {r for r in reactions_to_check
                                if abs(lp_model.variables.get(r).primal) - eps >= tolerance}

            if len(reactions_to_check & active_reactions) > 0:
                consistent_rxns |= active_reactions
                flipped = False
                singleton = False
            else:
                if singleton:
                    reactions_to_check -= rxns_subset
                if flipped:
                    flipped = False
                    singleton = True
                else:
                    flipped = True
                    if singleton:
                        if len(reactions_to_check) > 0:
                            rxn = sorted(reactions_to_check)[0]
                            rxns_subset = {rxn}
                        else:
                            rxns_subset = set()
                    else:
                        rxns_subset = reactions_to_check
                for rxn_id in rxns_subset:
                    reactions_dict[rxn_id] = {k: -v for k, v in
                                              reactions_dict[rxn_id].items()}
                    Fastcore.flip_variable(rxn_id, lp_model, reactions_dict[rxn_id])

        return consistent_rxns

    @staticmethod
    def fast_find_blocked(cobra_model, epsilon=1e-4, tolerance=1e-7, debug=False):
        consistent_reactions = Fastcore.fastcc(cobra_model, epsilon=epsilon,
                                               tolerance=tolerance, debug=debug)

        the_reactions = set([r.id for r in cobra_model.reactions])
        blocked_reactions = the_reactions - consistent_reactions
        return blocked_reactions

    # Private methods

    def _initialize_penalties(self, penalties):
        if penalties is None:
            penalties = {}
        for r in self._model.reactions:
            if r.id in self.core_reactions:
                self.penalties[r.id] = 0
            elif r.id in penalties:
                self.penalties[r.id] = penalties[r.id]
            else:
                self.penalties[r.id] = self._default_penalty

    def _consistency_check(self):
        # Finding the set of blocked reactions
        self.blocked_reactions = frozenset(Fastcore.fast_find_blocked(self._model))
        if len(self.blocked_reactions) > 0:
            self.gap_metabolites = []
            for m in self._model.metabolites:
                met_rxns_set = set([r.id for r in m.reactions])
                # If a metabolites participates only in blocked reaction
                # it is defined as a gap metabolite (Ponce-de-Leon et al. 2013)
                if len(met_rxns_set - self.blocked_reactions) == 0:
                    self.gap_metabolites.append(m.id)

    def _prune_inconsistencies(self):
        if len(self.blocked_reactions) > 0:
            blocked_reactions = [self._model.reactions.get_by_id(r) for r in self.blocked_reactions]
            self._model.remove_reactions(blocked_reactions)

        if len(self.gap_metabolites) > 0:
            gap_metabolites = [self._model.metabolites.get_by_id(m) for m in self.gap_metabolites]
            self._model.remove_metabolites(gap_metabolites)

    # Instance methods

    def fast_core(self):
        print("Running Fastcore")
        print("- Nº of core reactions %i " % len(self.consistent_core_reactions))
        reactions_dict = Fastcore.create_reactions_dict(self._model)
        consistent_subnet = set()
        self._consistent_sets = []

        all_rxns = set(self.consistent_reactions)
        core_rxns = set(self.consistent_core_reactions)
        non_core_rxns = all_rxns - core_rxns

        irreversible_rxns = {r for r in all_rxns if self.bounds[r][0] >= 0}
        irreversible_core_rxns = core_rxns & irreversible_rxns
        sparse_mode = self.find_sparse_mode(irreversible_core_rxns, non_core_rxns)

        if len(irreversible_core_rxns - sparse_mode):
            if self._debug_mode and len(irreversible_core_rxns - sparse_mode) > 0:
                print(irreversible_core_rxns - sparse_mode)
            raise Exception("Error: Inconsistent irreversible core reactions")

        consistent_subnet |= sparse_mode
        core_rxns -= consistent_subnet
        self._consistent_sets.append(consistent_subnet)
        if self._debug_mode:
            print("|Consistent Subnet|=", len(consistent_subnet))
            print("|Core Reactions|=", len(core_rxns))

        core_rxns_subset = {}
        flipped = False
        singleton = False
        while core_rxns:

            non_core_rxns -= consistent_subnet
            for r in consistent_subnet:
                self.penalties[r] = 0

            if singleton:
                rxn = sorted(core_rxns)[0]
                sparse_mode = self.find_sparse_mode({rxn}, non_core_rxns)
            else:
                sparse_mode = self.find_sparse_mode(core_rxns, non_core_rxns)

            consistent_subnet |= sparse_mode

            if self._debug_mode:
                print("|Consistent Subnet|=", len(consistent_subnet))

            if core_rxns & consistent_subnet:
                self._consistent_sets.append(consistent_subnet)
                core_rxns -= consistent_subnet
                flipped = False
                singleton = False
                if self._debug_mode:
                    print("|Core Reactions|=", len(core_rxns))
            else:
                if flipped:
                    flipped = False
                    singleton = True
                else:
                    flipped = True
                    if singleton:
                        rxn = sorted(core_rxns)[0]
                        core_rxns_subset = {rxn}
                    else:
                        core_rxns_subset = core_rxns

                for rxn_id in core_rxns_subset:
                    reactions_dict[rxn_id] = {k: -v for k, v in
                                              reactions_dict[rxn_id].items()}
                    Fastcore.flip_variable(rxn_id, self.LP7, reactions_dict[rxn_id])
                    Fastcore.flip_variable(rxn_id, self.LP10, reactions_dict[rxn_id])

        self._consistent_sets.append(consistent_subnet)
        print("Consistent sub-network found!")
        print("- Nº of reactions in consistent sub-network %i " % len(self._consistent_sets[-1]))

    def find_sparse_mode(self, core_reactions, non_core_reactions):

        assert len(core_reactions & non_core_reactions) == 0

        eps = 0.99 * self.epsilon

        if len(core_reactions) == 0:
            return set()

        Fastcore.update_optlang_lp7(self.LP7, core_reactions)
        self.LP7.optimize()
        if self.LP7.status == 'infeasible':
            raise Exception("INFEASIBLE LP7")

        active_core = {var.name for var in self.LP7.variables
                       if (var.name in core_reactions) and (var.primal - eps >= self.tolerance)}

        if len(active_core) == 0:
            return set()

        Fastcore.update_optlang_lp10(self.LP10, active_core, self.penalties, self.bounds,
                                     scaling_factor=self.scaling_factor, epsilon=self.epsilon)
        self.LP10.optimize()
        if self.LP10.status == 'infeasible':
            raise Exception("INFEASIBLE LP10")

        sparse_mode = []
        eps = 0.99 * self.epsilon
        for var in self.LP10.variables:
            if var.name.startswith(Fastcore.AUXILIARY_VAR_PREFIX):
                continue
            if abs(var.primal) - eps < self.tolerance:
                continue
            sparse_mode.append(var.name)

        return set(sparse_mode)

    def get_context_specific_set(self):
        if self._consistent_sets:
            return self._consistent_sets[-1]
        else:
            print("Before get a consistent network fast_core method should run")
            return None

    def build_context_specific_model(self):
        cs_rxns = self.get_context_specific_set()
        if cs_rxns:
            cs_model = self._original_model.copy()
            to_remove = []
            for r in cs_model.reactions:
                if r.id in cs_rxns:
                    continue
                to_remove.append(r)
            # Keep only reactions in the consistent subnetwork
            cs_model.remove_reactions(to_remove)
            # Removing metabolites not involved in any reaction
            cs_model.remove_metabolites([m for m in cs_model.metabolites
                                    if len(m.reactions) == 0])
        else:
            cs_model = None

        return cs_model



