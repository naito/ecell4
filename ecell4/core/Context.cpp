#include "Context.hpp"
#include <string>
#include <sstream>


namespace ecell4
{

std::pair<bool, MatchObject::context_type> uspmatch(
    const UnitSpecies& pttrn, const UnitSpecies& usp,
    const MatchObject::context_type& org)
{
    std::pair<bool, MatchObject::context_type>
        retval(std::make_pair(false, org));
    MatchObject::context_type& ctx(retval.second);

    if (rbex::is_wildcard(pttrn.name()))
    {
        if (rbex::is_pass_wildcard(pttrn.name()))
        {
            throw NotSupported(
                "A pass wildcard '_0' is not allowed to be a name of Species.");
        }
        else if (rbex::is_named_wildcard(pttrn.name()))
        {
            MatchObject::context_type::variable_container_type::const_iterator
                itr(ctx.globals.find(pttrn.name()));
            if (itr == ctx.globals.end())
            {
                ctx.globals[pttrn.name()] = usp.name();
            }
            else if ((*itr).second != usp.name())
            {
                return retval;
            }
        }
    }
    else if (pttrn.name() != usp.name())
    {
        return retval;
    }

    for (UnitSpecies::container_type::const_iterator j(pttrn.begin());
        j != pttrn.end(); ++j)
    {
        if (usp.has_site((*j).first))
        {
            const UnitSpecies::site_type& site(usp.get_site((*j).first));

            if ((*j).second.first != "")
            {
                if (site.first == "")
                {
                    return retval;
                }
                else if (rbex::is_pass_wildcard((*j).second.first))
                {
                    throw NotSupported(
                        "A pass wildcard '_0' is not allowed to be a state.");
                }
                else if (rbex::is_unnamed_wildcard((*j).second.first))
                {
                    ; // do nothing
                }
                else if (rbex::is_named_wildcard((*j).second.first))
                {
                    MatchObject::context_type::variable_container_type::const_iterator
                        itr(ctx.globals.find((*j).second.first));
                    if (itr == ctx.globals.end())
                    {
                        ctx.globals[(*j).second.first] = site.first;
                    }
                    else if ((*itr).second != site.first)
                    {
                        return retval;
                    }
                }
                else if ((*j).second.first != site.first)
                {
                    return retval;
                }
            }

            if (rbex::is_pass_wildcard((*j).second.second))
            {
                ; // just skip checking
            }
            else if ((*j).second.second == "")
            {
                if (site.second != "")
                {
                    return retval;
                }
            }
            else
            {
                if (site.second == "")
                {
                    return retval;
                }
                else if (rbex::is_unnamed_wildcard((*j).second.second))
                {
                    continue;
                }
                else if (rbex::is_named_wildcard((*j).second.second))
                {
                    throw NotSupported(
                        "A named wildcard is not allowed to be a bond.");
                }

                MatchObject::context_type::variable_container_type::const_iterator
                    itr(ctx.locals.find((*j).second.second));
                if (itr == ctx.locals.end())
                {
                    ctx.locals[(*j).second.second] = site.second;
                }
                else if ((*itr).second != site.second)
                {
                    return retval;
                }

            }
        }
        else
        {
            return retval;
        }
    }

    retval.first = true;
    return retval;
}

bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx)
{
    if (itr == end)
    {
        // for (MatchObject::context_type::iterator_container_type::const_iterator
        //     i(ctx.iterators.begin()); i != ctx.iterators.end(); ++i)
        //     std::cout << *i << " ";
        // std::cout << std::endl;
        return true;
    }

    MatchObject obj(*itr);
    ++itr;

    std::pair<bool, MatchObject::context_type> retval(obj.match(sp, ctx));
    while (retval.first)
    {
        if (__spmatch(itr, end, sp, retval.second))
        {
            return true;
        }
        retval = obj.next();
    }
    return false;
}

bool spmatch(const Species& pttrn, const Species& sp)
{
    SpeciesExpressionMatcher sexp(pttrn);
    return sexp.match(sp);
    // MatchObject::context_type ctx;
    // return __spmatch(pttrn.begin(), pttrn.end(), sp, ctx);
}

Integer count_spmatches(const Species& pttrn, const Species& sp)
{
    MatchObject::context_type::variable_container_type globals;
    return count_spmatches(pttrn, sp, globals);
}

Integer count_spmatches(const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals)
{
    SpeciesExpressionMatcher sexp(pttrn);
    if (!sexp.match(sp, globals))
    {
        return 0;
    }
    Integer n(1);
    while (sexp.next())
    {
        ++n;
    }
    return n;
}

Integer count_spmatches(const std::vector<UnitSpecies>& pttrn, const std::vector<UnitSpecies>& sp)
{
    MatchObject::context_type::variable_container_type globals;
    return count_spmatches(pttrn, sp, globals);
}

Integer count_spmatches(const std::vector<UnitSpecies>& pttrn, const std::vector<UnitSpecies>& sp,
    const MatchObject::context_type::variable_container_type& globals)
{
    SpeciesExpressionMatcher sexp(pttrn);
    if (!sexp.match(sp, globals))
    {
        return 0;
    }
    Integer n(1);
    while (sexp.next())
    {
        ++n;
    }
    return n;
}

std::pair<bool, MatchObject::context_type> __rrmatch(
    const ReactionRule& rr,
    const ReactionRule::reactant_container_type& reactants,
    const MatchObject::context_type::variable_container_type& globals,
    ReactionRule::reactant_container_type::const_iterator i,
    ReactionRule::reactant_container_type::const_iterator j)
{
    SpeciesExpressionMatcher m(*i);
    if (!m.match(*j, globals))
    {
        return std::make_pair(false, MatchObject::context_type());
    }

    ++i;
    ++j;
    if (i == rr.reactants().end() || j == reactants.end())
    {
        return std::make_pair(true, m.context());
    }

    do
    {
        if (__rrmatch(rr, reactants, m.context().globals, i, j).first)
        {
            return std::make_pair(true, m.context());
        }
    } while (m.next());
    return std::make_pair(false, MatchObject::context_type());
}

std::pair<bool, MatchObject::context_type> MatchObject::next()
{
    std::vector<UnitSpecies>::const_iterator itr_start = target_.begin();
    for (; itr_ != target_.end(); ++itr_)
    {
        const Species::container_type::difference_type
            pos(distance(itr_start, itr_));
        if (std::find(ctx_.iterators.begin(), ctx_.iterators.end(), pos)
            != ctx_.iterators.end())
        {
            continue;
        }

        const UnitSpecies& usp(*itr_);
        std::pair<bool, MatchObject::context_type>
            retval(uspmatch(pttrn_, usp, ctx_));
        if (retval.first)
        {
            retval.second.iterators.push_back(pos);
            ++itr_;
            return retval;
        }
    }
    return std::make_pair(false, MatchObject::context_type());
}

unsigned int concatenate_units(std::vector<UnitSpecies>& units1, const std::vector<UnitSpecies>& units2, const unsigned int bond_stride)
{
    // const std::vector<UnitSpecies> units2 = sp.units();
    units1.reserve(units1.size() + units2.size());

    utils::get_mapper_mf<std::string, std::string>::type bond_cache;

    for (std::vector<UnitSpecies>::const_iterator j(units2.begin());
        j != units2.end(); ++j)
    {
        units1.push_back(*j);

        for (UnitSpecies::container_type::const_iterator
            k((*j).begin()); k != (*j).end(); ++k)
        {
            const std::string& bond((*k).second.second);
            if (bond != "" && !rbex::is_wildcard(bond))
            {
                utils::get_mapper_mf<std::string, std::string>::type::const_iterator
                    it = bond_cache.find(bond);
                if (it == bond_cache.end())
                {
                    const std::string newbond = itos(bond_stride + 1 + bond_cache.size());
                    bond_cache.insert(std::make_pair(bond, newbond));
                    units1.back().at(std::distance((*j).begin(), k)).second.second = newbond;
                }
                else
                {
                    units1.back().at(std::distance((*j).begin(), k)).second.second = (*it).second;
                }
            }
        }
    }
    return bond_cache.size();
}

bool is_correspondent(const UnitSpecies& usp1, const UnitSpecies& usp2)
{
    if (usp1.name() != usp2.name() || usp1.num_sites() != usp2.num_sites())
    {
        return false;
    }

    for (UnitSpecies::container_type::const_iterator
        i(usp1.begin()), j(usp2.begin()); i != usp1.end() && j != usp2.end(); ++i, ++j)
    {
        if ((*i).first != (*j).first)
        {
            return false;
        }
    }
    return true;
}

unsigned int tag_units(
    std::vector<unsigned int>& groups,
    const unsigned int& group_id,
    const unsigned int& idx,
    const std::vector<std::vector<std::vector<UnitSpecies>::size_type> >& connections,
    const unsigned int& notyet)
{
    if (groups[idx] != notyet)
    {
        // assert(groups[idx] < group_id);
        return group_id;
    }

    groups[idx] = group_id;

    for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
        i(connections[idx].begin()); i != connections[idx].end(); ++i)
    {
        const unsigned int _gid = tag_units(groups, group_id, *i, connections, notyet);
        // assert(_gid == group_id);
    }

    return group_id + 1;
}

DummyReactionRule format_reaction_rule(const DummyReactionRule& rr)
{
    std::vector<Dummy> reactants;
    reactants.reserve(rr.reactants().size());
    for (std::vector<Dummy>::const_iterator i(rr.reactants().begin());
        i != rr.reactants().end(); ++i)
    {
        reactants.push_back(Dummy(format_species((*i).units())));
    }

    std::vector<Dummy> products;
    products.reserve(rr.products().size());
    for (std::vector<Dummy>::const_iterator i(rr.products().begin());
        i != rr.products().end(); ++i)
    {
        products.push_back(Dummy(format_species((*i).units())));
    }

    std::sort(reactants.begin(), reactants.end());
    std::sort(products.begin(), products.end());

    return rr(reactants, products);
}

namespace rbex_extras
{

std::vector<DummyReactionRule> generate_reaction_rules(
    const DummyReactionRule& org, const Dummy& sp1)
{
    std::vector<Dummy> reactants(1, sp1);
    return org.generate(reactants);
}

std::vector<DummyReactionRule> generate_reaction_rules(
    const DummyReactionRule& org, const Dummy& sp1, const Dummy& sp2)
{
    if (org.reactants().size() != 2)
    {
        return std::vector<DummyReactionRule>(0);
    }

    std::vector<Dummy> reactants(2);
    reactants[0] = sp1;
    reactants[1] = sp2;
    std::vector<DummyReactionRule> res = org.generate(reactants);

    // if (org.reactants()[0] != org.reactants()[1])
    // {
    //     reactants[0] = sp2;
    //     reactants[1] = sp1;
    //     std::vector<DummyReactionRule> _res = org.generate(reactants);
    //     std::copy(_res.begin(), _res.end(), back_inserter(res));
    // }
    return res;
}

bool check_stoichiometry(const Dummy& sp,
    const std::map<Dummy, Integer>& max_stoich)
{
    for (std::map<Dummy, Integer>::const_iterator i(max_stoich.begin());
        i != max_stoich.end(); ++i)
    {
        if (count_spmatches((*i).first.units(), sp.units()) > (*i).second)
        {
            return false;
        }
    }
    return true;
}

bool check_stoichiometry(const DummyReactionRule& rr,
    const std::map<Dummy, Integer>& max_stoich)
{
    for (std::vector<Dummy>::const_iterator
        i(rr.products().begin()); i != rr.products().end(); ++i)
    {
        if (!check_stoichiometry(*i, max_stoich))
        {
            return false;
        }
    }
    return true;
}

void __add_reaction_rules(
    const std::vector<DummyReactionRule>& reaction_rules,
    std::vector<DummyReactionRule>& reactions, std::vector<Dummy>& newseeds,
    const std::vector<Dummy>& seeds,
    const std::map<Dummy, Integer>& max_stoich)
{
    for (std::vector<DummyReactionRule>::const_iterator i(reaction_rules.begin());
        i != reaction_rules.end(); ++i)
    {
        const DummyReactionRule& rr(*i);
        if (!check_stoichiometry(rr, max_stoich))
        {
            continue;
        }

        reactions.push_back(rr);

        for (std::vector<Dummy>::const_iterator
            j(rr.products().begin()); j != rr.products().end(); ++j)
        {
            const Dummy sp(format_species((*j).units()));
            if (std::find(newseeds.begin(), newseeds.end(), sp)
                == newseeds.end()
                && std::find(seeds.begin(), seeds.end(), sp)
                == seeds.end())
            {
                newseeds.push_back(sp);
            }
        }
    }
}

void __apply_reaction_rules(
    std::vector<Dummy>& seeds,
    const std::vector<DummyReactionRule>& rules,
    std::vector<DummyReactionRule>& reactions,
    std::vector<Dummy>& allseeds,
    const std::map<Dummy, Integer>& max_stoich,
    std::vector<std::vector<std::vector<Dummy>::size_type> >& candidates)
{
    std::vector<Dummy>::size_type stride = allseeds.size();

    std::vector<std::vector<std::vector<Dummy>::size_type> >::size_type count = 0;
    for (std::vector<DummyReactionRule>::const_iterator
        i(rules.begin()); i != rules.end(); ++i)
    {
        const DummyReactionRule& rr(*i);

        switch (rr.reactants().size())
        {
        case 0:
            continue;
        case 1:
            continue;
        case 2:
            {
                if (candidates.size() < count + 2 + 1)
                {
                    candidates.resize(count + 2 + 1);
                }

                SpeciesExpressionMatcher reactant0(rr.reactants()[0].units());
                SpeciesExpressionMatcher reactant1(rr.reactants()[1].units());

                std::vector<Dummy>::size_type idx = stride;
                for (std::vector<Dummy>::const_iterator j(seeds.begin());
                    j != seeds.end(); ++j, ++idx)
                {
                    if (reactant0.match((*j).units()))
                    {
                        candidates[count].push_back(idx);
                    }
                    if (reactant1.match((*j).units()))
                    {
                        candidates[count + 1].push_back(idx);
                    }
                }
                count += 2;
            }
            break;
        default:
            throw NotImplemented(
                "No reaction rule with more than two reactants is accepted.");
        }
    }

    // allseeds.insert(allseeds.begin(), seeds.begin(), seeds.end());
    allseeds.insert(allseeds.end(), seeds.begin(), seeds.end());

    std::vector<Dummy> newseeds;

    // std::cout << "Start Generations: " << stride << ", " << seeds.size() << std::endl;

    count = 0;
    for (std::vector<DummyReactionRule>::const_iterator
        i(rules.begin()); i != rules.end(); ++i)
    {
        const DummyReactionRule& rr(*i);

        switch (rr.reactants().size())
        {
        case 0:
            continue;
        case 1:
            for (std::vector<Dummy>::const_iterator j(seeds.begin());
                j != seeds.end(); ++j)
            {
                __add_reaction_rules(
                    generate_reaction_rules(rr, *j),
                    reactions, newseeds, allseeds, max_stoich);
            }
            break;
        case 2:
            {
                for (std::vector<std::vector<Dummy>::size_type>::const_iterator j(candidates[count].begin());
                    j != candidates[count].end(); ++j)
                {
                    for (std::vector<std::vector<Dummy>::size_type>::const_iterator k(candidates[count + 1].begin());
                        k != candidates[count + 1].end(); ++k)
                    {
                        // std::cout << "Reaction [" << count << "]: " << (*j) << ", " << (*k) << std::endl;

                        if ((*j) < stride && (*k) < stride)
                        {
                            continue;
                        }

                        // std::cout << "Generate Reaction [" << count << ":" << rr.data.as_string() << "]: ";

                        const Dummy& reactant0(allseeds[(*j)]);
                        const Dummy& reactant1(allseeds[(*k)]);

                        // std::cout << reactant0.serial() << " + " << reactant1.serial() << std::endl;

                        __add_reaction_rules(
                            generate_reaction_rules(rr, reactant0, reactant1),
                            reactions, newseeds, allseeds, max_stoich);
                    }
                }

                // for (std::vector<Dummy>::const_iterator j(seeds.begin());
                //     j != seeds.end(); ++j)
                // {
                //     const std::vector<Dummy>::const_iterator start(
                //         allseeds.begin()
                //         + std::distance<std::vector<Dummy>::const_iterator>(
                //             seeds.begin(), j));
                //     for (std::vector<Dummy>::const_iterator
                //         k(start); k != allseeds.end(); ++k)
                //     {
                //         __add_reaction_rules(
                //             generate_reaction_rules(rr, *j, *k),
                //             reactions, newseeds, allseeds, max_stoich);
                //     }
                // }

                count += 2;
            }
            break;
        default:
            throw NotImplemented(
                "No reaction rule with more than two reactants is accepted.");
        }
    }

    // for (std::vector<Dummy>::const_iterator j(newseeds.begin());
    //     j != newseeds.end(); ++j)
    // {
    //     std::cout << "NEW: " << (*j).serial() << std::endl;
    // }
    seeds.swap(newseeds);
}

std::pair<std::vector<ReactionRule>, bool> apply_reaction_rules(
    std::vector<Species>& allseeds,
    const std::vector<ReactionRule>& rules,
    const Integer max_itr,
    const std::map<Species, Integer>& max_stoich)
{
    std::vector<Dummy> _allseeds(0);
    std::vector<Dummy> _seeds = DummyBinder<Species, Dummy>::encode(allseeds);
    std::vector<DummyReactionRule> _rules = DummyBinder<ReactionRule, DummyReactionRule>::encode(rules);
    std::vector<DummyReactionRule> _reactions;

    std::map<Dummy, Integer> _max_stoich;
    for (std::map<Species, Integer>::const_iterator i(max_stoich.begin());
        i != max_stoich.end(); ++i)
    {
        _max_stoich[DummyBinder<Species, Dummy>::encode((*i).first)] = (*i).second;
    }

    for (std::vector<DummyReactionRule>::const_iterator
        i(_rules.begin()); i != _rules.end(); ++i)
    {
        const DummyReactionRule& rr(*i);
        if (rr.reactants().size() == 0 && check_stoichiometry(rr, _max_stoich))
        {
            _reactions.push_back(rr);
            for (std::vector<Dummy>::const_iterator
                j(rr.products().begin()); j != rr.products().end(); ++j)
            {
                const Dummy sp(format_species((*j).units()));
                if (std::find(_seeds.begin(), _seeds.end(), sp)
                    == _seeds.end())
                {
                    _seeds.push_back(sp);
                }
            }
        }
    }

    {
        std::vector<std::vector<std::vector<Dummy>::size_type> > _candidates;

        Integer iteration_count(0);
        while (_seeds.size() > 0 && iteration_count < max_itr)
        {
            __apply_reaction_rules(_seeds, _rules, _reactions, _allseeds, _max_stoich, _candidates);
            iteration_count += 1;
        }
    }

    const bool is_completed = (_seeds.size() == 0);
    if (!is_completed)
    {
        _allseeds.insert(_allseeds.begin(), _seeds.begin(), _seeds.end());
    }

    std::transform(_reactions.begin(), _reactions.end(), _reactions.begin(), (DummyReactionRule (*)(const DummyReactionRule&))format_reaction_rule);

    allseeds = DummyBinder<Species, Dummy>::decode(_allseeds);
    return std::make_pair(DummyBinder<ReactionRule, DummyReactionRule>::decode(_reactions), is_completed);
}

} // rbm

} // ecell4
