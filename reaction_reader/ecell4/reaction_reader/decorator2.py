import species
import copy
import decorator
import parseobj
import functools


def generate_Species(obj):
    if isinstance(obj, parseobj.AnyCallable):
        obj = obj._as_ParseObj()

    if isinstance(obj, parseobj.ParseObj):
        sp = species.Species()
        for elem in obj._elements():
            su = species.Subunit(elem.name)
            if elem.args is not None:
                for mod in elem.args:
                    if (not (isinstance(mod, parseobj.ParseObj)
                            or isinstance(mod, parseobj.AnyCallable))
                        or mod._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (mod))
                    arg = mod._elements()[0]
                    name, binding = arg.name, arg.modification
                    if binding is None:
                        su.add_modification(name, "", "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, "", binding)
            if elem.kwargs is not None:
                for name, value in elem.kwargs.items():
                    if (not (isinstance(value, parseobj.ParseObj)
                            or isinstance(value, parseobj.AnyCallable))
                        or value._size() != 1):
                        raise RuntimeError, (
                            "invalid argument [%s] found." % (mod))
                    arg = value._elements()[0]
                    state, binding = str(arg.name), arg.modification
                    if binding is None:
                        su.add_modification(name, state, "")
                    else:
                        binding = str(binding)
                        if not (binding.isdigit() or binding == ""
                            or binding[0] == "_"):
                            raise RuntimeError, (
                                "invalid binding [%s] given." % (binding))
                        su.add_modification(name, state, binding)
            sp.add_subunit(su)
        return (sp, )
    elif isinstance(obj, parseobj.InvExp):
        return (None, )
    elif isinstance(obj, parseobj.AddExp):
        subobjs = obj._elements()
        return tuple(generate_Species(subobj)[0] for subobj in subobjs)

    raise RuntimeError, 'invalid expression; "%s" given' % str(obj)

def generate_ReactionRule(lhs, rhs, k=0.0):
    # if len(lhs) == 0:
    #     if len(rhs) != 1:
    #         raise RuntimeError, (
    #             "the number of products must be 1; %d given" % len(rhs))
    #     return ecell4.core.create_synthesis_reaction_rule(rhs[0], k)
    if len(lhs) == 1 or len(lhs) == 2:
        return species.ReactionRule(lhs, rhs)
    raise RuntimeError, (
        "the number of reactants must be less than 3; %d given" % len(lhs))

class SpeciesAttributesCallback(decorator.Callback):

    def __init__(self, *args):
        decorator.Callback.__init__(self)

        self.bitwise_operations = []

    def get(self):
        return copy.copy(self.bitwise_operations)

    def notify_bitwise_operations(self, obj):
        if not isinstance(obj, parseobj.OrExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))

        lhs, rhs = obj._lhs, obj._rhs

        species_list = generate_Species(lhs)
        if len(species_list) != 1:
            raise RuntimeError, (
                'only a single species must be given; %d given'
                % len(species_list))

        sp = species_list[0]
        if sp is None:
            raise RuntimeError, "no species given [%s]" % (repr(obj))

        self.bitwise_operations.append((sp, rhs))

    def notify_comparisons(self, obj):
        raise RuntimeError, (
            'ReactionRule definitions are not allowed'
            + ' in "species_attributes"')

class ReactionRulesCallback(decorator.Callback):

    def __init__(self):
        decorator.Callback.__init__(self)

        self.comparisons = []

    def get(self):
        return copy.copy(self.comparisons)

    def notify_comparisons(self, obj):
        if not isinstance(obj, parseobj.CmpExp):
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))
        elif isinstance(obj, parseobj.NeExp):
            warnings.warn('"<>" is deprecated; use "==" instead',
                          DeprecationWarning)

        lhs, rhs = obj._lhs, obj._rhs

        if isinstance(lhs, parseobj.OrExp):
            lhs = lhs._lhs

        if isinstance(rhs, parseobj.OrExp):
            rhs = rhs._lhs

        lhs, rhs = generate_Species(lhs), generate_Species(rhs)
        lhs = tuple(sp for sp in lhs if sp is not None)
        rhs = tuple(sp for sp in rhs if sp is not None)

        if isinstance(obj, parseobj.EqExp) or isinstance(obj, parseobj.NeExp):
            self.comparisons.append(generate_ReactionRule(lhs, rhs))
            self.comparisons.append(generate_ReactionRule(rhs, lhs))
        elif isinstance(obj, parseobj.GtExp):
            self.comparisons.append(generate_ReactionRule(lhs, rhs))
        else:
            raise RuntimeError, 'an invalid object was given [%s]' % (repr(obj))

species_attributes = functools.partial(decorator.parse_decorator, SpeciesAttributesCallback)
reaction_rules = functools.partial(decorator.parse_decorator, ReactionRulesCallback)