#ifndef ECELL4_CONTEXT_HPP
#define ECELL4_CONTEXT_HPP

#include "get_mapper_mf.hpp"
#include "Species.hpp"
#include "ReactionRule.hpp"
#include <boost/array.hpp>


namespace ecell4
{

namespace rbex
{

inline bool is_wildcard(const std::string& name)
{
    return (name.size() > 0 && name[0] == '_');
}

inline bool is_unnamed_wildcard(const std::string& name)
{
    return name == "_";
}

inline bool is_pass_wildcard(const std::string& name)
{
    return name == "_0";
}

inline bool is_named_wildcard(const std::string& name)
{
    return (name.size() > 1 && name[0] == '_' && !is_pass_wildcard(name));
}

} // rbex

class MatchObject
{
public:

    typedef struct
    {
        typedef std::vector<Species::container_type::difference_type>
            iterator_container_type;
        typedef utils::get_mapper_mf<std::string, std::string>::type
            variable_container_type;

        iterator_container_type iterators;
        variable_container_type locals;
        variable_container_type globals;
    } context_type;

public:

    MatchObject(const UnitSpecies& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    virtual ~MatchObject()
    {
        ;
    }

    std::pair<bool, context_type> match(
        const Species& sp, const context_type& ctx)
    {
        // target_ = sp;
        target_ = sp.units();
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> match(
        const std::vector<UnitSpecies>& target, const context_type& ctx)
    {
        target_ = target;
        itr_ = target_.begin();
        ctx_ = ctx;
        return next();
    }

    std::pair<bool, context_type> next();

protected:

    UnitSpecies pttrn_;
    // Species target_;
    std::vector<UnitSpecies> target_;
    Species::container_type::const_iterator itr_;
    context_type ctx_;
};

std::pair<bool, MatchObject::context_type>
uspmatch(const UnitSpecies& pttrn, const UnitSpecies& sp,
    const MatchObject::context_type& org);
bool __spmatch(
    Species::container_type::const_iterator itr,
    const Species::container_type::const_iterator& end,
    const Species& sp, const MatchObject::context_type& ctx);
bool spmatch(const Species& pttrn, const Species& sp);
Integer count_spmatches(const Species& pttrn, const Species& sp);
Integer count_spmatches(
    const Species& pttrn, const Species& sp,
    const MatchObject::context_type::variable_container_type& globals);
Integer count_spmatches(const std::vector<UnitSpecies>& pttrn, const std::vector<UnitSpecies>& sp);
Integer count_spmatches(
    const std::vector<UnitSpecies>& pttrn, const std::vector<UnitSpecies>& sp,
    const MatchObject::context_type::variable_container_type& globals);

// ReactionRule create_reaction_rule_formatted(
//     const ReactionRule::reactant_container_type& reactants,
//     const ReactionRule::product_container_type& products, const Real k);

// inline ReactionRule create_reaction_rule_formatted(const ReactionRule& rr)
// {
//     return create_reaction_rule_formatted(rr.reactants(), rr.products(), rr.k());
// }

class SpeciesExpressionMatcher
{
public:

    typedef MatchObject::context_type context_type;

public:

    SpeciesExpressionMatcher(const Species& pttrn)
        : pttrn_(pttrn.units())
    {
        ;
    }

    SpeciesExpressionMatcher(const std::vector<UnitSpecies>& pttrn)
        : pttrn_(pttrn)
    {
        ;
    }

    virtual ~SpeciesExpressionMatcher()
    {
        ;
    }

    bool match(const Species& sp)
    {
        context_type::variable_container_type globals;
        return match(sp, globals);
    }

    bool match(
        const Species& sp, const context_type::variable_container_type& globals)
    {
        matches_.clear();
        for (Species::container_type::const_iterator i(pttrn_.begin());
            i != pttrn_.end(); ++i)
        {
            matches_.push_back(MatchObject(*i));
        }

        // target_ = sp;
        target_ = sp.units();
        itr_ = matches_.begin();
        context_type ctx;
        ctx.globals = globals;
        return __match(ctx);
    }

    bool match(const std::vector<UnitSpecies>& sp)
    {
        context_type::variable_container_type globals;
        return match(sp, globals);
    }

    bool match(
        const std::vector<UnitSpecies>& sp, const context_type::variable_container_type& globals)
    {
        matches_.clear();
        for (Species::container_type::const_iterator i(pttrn_.begin());
            i != pttrn_.end(); ++i)
        {
            matches_.push_back(MatchObject(*i));
        }

        target_ = sp;
        itr_ = matches_.begin();
        context_type ctx;
        ctx.globals = globals;
        return __match(ctx);
    }

    bool __match(const context_type& ctx)
    {
        if (itr_ == matches_.end())
        {
            ctx_ = ctx;
            return true;
        }

        std::pair<bool, context_type> retval((*itr_).match(target_, ctx));
        while (retval.first)
        {
            ++itr_;
            const bool succeeded(__match(retval.second));
            if (succeeded)
            {
                return true;
            }
            --itr_;
            retval = (*itr_).next();
        }
        return false;
    }

    bool next()
    {
        if (itr_ != matches_.end())
        {
            return false;
        }
        else if (matches_.size() == 0)
        {
            return true;
        }

        do
        {
            --itr_;
            std::pair<bool, context_type> retval((*itr_).next());
            while (retval.first)
            {
                ++itr_;
                const bool succeeded(__match(retval.second));
                if (succeeded)
                {
                    return true;
                }
                --itr_;
                retval = (*itr_).next();
            }
        }
        while (itr_ != matches_.begin());
        return false;
    }

    Integer count(const Species& sp)
    {
        context_type::variable_container_type globals;
        if (!match(sp, globals))
        {
            return 0;
        }
        Integer n(1);
        while (next())
        {
            ++n;
        }
        return n;
    }

    const context_type& context() const
    {
        return ctx_;
    }

protected:

    // Species pttrn_;
    // Species target_;
    std::vector<UnitSpecies> pttrn_;
    std::vector<UnitSpecies> target_;
    std::vector<MatchObject> matches_;
    std::vector<MatchObject>::iterator itr_;
    context_type ctx_;
};

inline std::string itos(unsigned int val)
{
    std::stringstream ss;
    ss << val;
    return ss.str();
}

unsigned int concatenate_units(std::vector<UnitSpecies>& units1, const std::vector<UnitSpecies>& units2, const unsigned int bond_stride);
bool is_correspondent(const UnitSpecies& usp1, const UnitSpecies& usp2);
unsigned int tag_units(
    std::vector<unsigned int>& groups,
    const unsigned int& group_id,
    const unsigned int& idx,
    const std::vector<std::vector<std::vector<UnitSpecies>::size_type> >& connections,
    const unsigned int& notyet);

template <typename Element>
std::vector<Element> group_units(
    const std::vector<UnitSpecies>& units,
    const std::vector<unsigned int>& groups, const unsigned int num_groups)
{
    // 8. Divide units into Species

    std::vector<Element> products;
    products.resize(num_groups);

    // {
    //     std::vector<unsigned int>::size_type idx = 0;
    //     for (std::vector<UnitSpecies>::iterator i(units.begin());
    //         i != units.end(); ++i; ++idx)
    //     {
    //         products[groups[idx]].add_unit(*i);
    //     }
    // }

    for (unsigned int idx(0); idx != num_groups; ++idx)
    {
        utils::get_mapper_mf<std::string, std::string>::type new_bonds;
        unsigned int stride(1);

        for (std::vector<unsigned int>::const_iterator
            i(groups.begin()); i != groups.end(); ++i)
        {
            if (idx != *i)
            {
                continue;
            }

            UnitSpecies usp(units[std::distance(groups.begin(), i)]); //XXX: copy

            if (usp == UnitSpecies())
            {
                continue;
            }

            for (UnitSpecies::container_type::size_type j(0);
                j != usp.num_sites(); ++j)
            {
                UnitSpecies::container_type::value_type&
                    site(usp.at(j));
                const std::string bond(site.second.second);
                if (bond == "" || rbex::is_wildcard(bond))
                {
                    continue;
                }

                utils::get_mapper_mf<std::string, std::string>::type::const_iterator
                    itr(new_bonds.find(bond));
                if (itr == new_bonds.end())
                {
                    const std::string new_bond(itos(stride));
                    ++stride;
                    new_bonds[bond] = new_bond;
                    site.second.second = new_bond;
                }
                else
                {
                    site.second.second = (*itr).second;
                }
            }

            products[idx].add_unit(usp);
        }

        // products[idx] = format_species(products[idx]);
    }

    products.erase(std::remove(products.begin(), products.end(), Element()), products.end());

    return products;
}

template <typename First, typename Second>
struct DummyBinder
{
    static std::vector<Second> encode(std::vector<First> const& elements)
    {
        std::vector<Second> _elements;
        _elements.reserve(elements.size());
        for (typename std::vector<First>::const_iterator i(elements.begin());
            i != elements.end(); ++i)
        {
            _elements.push_back(encode(*i));
        }
        return _elements;
    }

    static std::vector<First> decode(std::vector<Second> const& elements)
    {
        std::vector<First> _elements;
        _elements.reserve(elements.size());
        for (typename std::vector<Second>::const_iterator i(elements.begin());
            i != elements.end(); ++i)
        {
            _elements.push_back(decode(*i));
        }
        return _elements;
    }

    static Second encode(First const& elem)
    {
        return Second(elem);
    }

    static First decode(Second const& elem)
    {
        return elem.get();
    }
};

template <typename Element>
class ReactionRuleExpressionMatcher
{
public:

    typedef Element element_type;

    typedef MatchObject::context_type context_type;
    // typedef ReactionRule::reactant_container_type reactant_container_type;
    typedef std::vector<element_type> element_container_type;

public:

    ReactionRuleExpressionMatcher(const element_container_type& reactants, const element_container_type& products, const ReactionRule::policy_type& policy)
        : pttrn_reactants_(reactants), pttrn_products_(products), policy_(policy)
    {
        ;
    }

    virtual ~ReactionRuleExpressionMatcher()
    {
        ;
    }

    bool match(const element_type& sp)
    {
        element_container_type reactants;
        reactants.push_back(sp);
        // permutation_.clear();
        // permutation_.push_back(0);
        return __match(reactants);
    }

    bool match(const element_type& sp1, const element_type& sp2)
    {
        element_container_type reactants;
        reactants.push_back(sp1);
        reactants.push_back(sp2);
        // permutation_.clear();
        // permutation_.push_back(0);
        // permutation_.push_back(1);
        return __match(reactants);
    }

    // bool match(const element_container_type& reactants,
    //            const std::vector<typename element_container_type::size_type>& permutation)
    // {
    //     permutation_ = permutation;
    //     return __match(reactants);
    // }

    bool match(const element_container_type& reactants)
    {
        // permutation_.clear();
        // permutation_.reserve(reactants.size());
        // for (typename std::vector<typename element_container_type::size_type>::size_type i(0);
        //     i != reactants.size(); ++i)
        // {
        //     permutation_.push_back(i);
        // }
        return __match(reactants);
    }

    bool __match(const element_container_type& reactants)
    {
        if (pttrn_reactants_.size() != reactants.size())
        {
            return false;
        }

        matchers_.clear();
        for (typename element_container_type::const_iterator
            i(pttrn_reactants_.begin()); i != pttrn_reactants_.end(); ++i)
        {
            matchers_.push_back(SpeciesExpressionMatcher((*i).units()));
        }

        target_ = reactants; //XXX: copy?
        itr_ = matchers_.begin();
        context_type::variable_container_type globals;
        return __submatch(globals);
    }

    bool __submatch(const context_type::variable_container_type& globals)
    {
        if (itr_ == matchers_.end())
        {
            return true;
        }

        bool retval((*itr_).match(
            target_[std::distance(matchers_.begin(), itr_)].units(),
            globals));
        while (retval)
        {
            const context_type::variable_container_type&
                globals_prev((*itr_).context().globals);
            ++itr_;
            const bool succeeded(__submatch(globals_prev));
            if (succeeded)
            {
                return true;
            }
            --itr_;
            retval = (*itr_).next();
        }
        return false;
    }

    bool next()
    {
        if (itr_ != matchers_.end() || pttrn_reactants_.size() == 0)
        {
            return false;
        }
        else if (matchers_.size() == 0)
        {
            return true;
        }

        do
        {
            --itr_;
            bool retval((*itr_).next());
            while (retval)
            {
                const context_type::variable_container_type&
                    globals_prev((*itr_).context().globals);
                ++itr_;
                const bool succeeded(__submatch(globals_prev));
                if (succeeded)
                {
                    return true;
                }
                --itr_;
                retval = (*itr_).next();
            }
        }
        while (itr_ != matchers_.begin());
        return false;
    }

    std::pair<bool, context_type> __match(
        const context_type::variable_container_type& globals,
        typename element_container_type::const_iterator i,
        typename element_container_type::const_iterator j)
    {
        SpeciesExpressionMatcher m((*i).units());
        if (!m.match(*j, globals))
        {
            return std::make_pair(false, context_type());
        }

        ++i;
        ++j;
        if (i == pttrn_reactants_.end() || j == target_.end())
        {
            return std::make_pair(true, m.context());
        }

        do
        {
            if (__match(m.context().globals, i, j).first)
            {
                return std::make_pair(true, m.context());
            }
        } while (m.next());
        return std::make_pair(false, context_type());
    }

    context_type context() const
    {
        context_type ctx;
        if (matchers_.size() == 0)
        {
            return ctx;
        }

        ctx.globals = matchers_.back().context().globals;

        std::vector<unsigned int> strides;
        strides.reserve(target_.size());
        {
            unsigned int stride = 0;
            for (typename element_container_type::const_iterator
                i(target_.begin()); i != target_.end(); ++i)
            {
                strides.push_back(stride);
                stride += (*i).units().size();
            }
        }

        for (std::vector<SpeciesExpressionMatcher>::const_iterator
            i(matchers_.begin()); i != matchers_.end(); ++i)
        {
            const unsigned int idx1 = std::distance(matchers_.begin(), i);  // a position in matcher_
            // const unsigned int idx2 = permutation_[idx1];  // a position in reactants
            const unsigned int idx2 = idx1;  // a position in reactants
            const unsigned int stride = strides[idx2];

            for (context_type::iterator_container_type::const_iterator
                j((*i).context().iterators.begin());
                j != (*i).context().iterators.end(); ++j)
            {
                // const unsigned int idx3 = std::distance((*i).context().iterators.begin(), j);  // a position in context.iterators
                const unsigned int idx4 = (*j);  // a position in units of a Species

                ctx.iterators.push_back(idx4 + stride);
            }
        }

        return ctx;
    }

    typedef struct
    {
        std::vector<UnitSpecies> products;
        std::vector<std::vector<UnitSpecies>::size_type> correspo;
        std::vector<std::vector<UnitSpecies>::size_type> removed;
        std::vector<UnitSpecies>::size_type reserved;
    } operation_type;

    operation_type compile()
    {
        typedef std::vector<UnitSpecies>::size_type size_type;
        typedef std::vector<UnitSpecies>::const_iterator const_iterator;

        operation_type res = {};
        std::vector<UnitSpecies>& products = res.products;
        std::vector<size_type>& correspo = res.correspo;
        std::vector<size_type>& removed = res.removed;

        // 1. Concatenate units of a pattern

        std::vector<UnitSpecies> reactants;
        for (typename element_container_type::const_iterator
            i(pttrn_reactants_.begin()); i != pttrn_reactants_.end(); ++i)
        {
            std::vector<UnitSpecies> const units = (*i).units();
            reactants.reserve(reactants.size() + units.size());
            std::copy(units.begin(), units.end(), std::back_inserter(reactants));
        }

        res.reserved = reactants.size();

        int product_bond_stride = 0;
        for (typename element_container_type::const_iterator
            i(pttrn_products_.begin()); i != pttrn_products_.end(); ++i)
        {
            product_bond_stride += concatenate_units(products, (*i).units(), product_bond_stride);
        }

        // 2. Check correspondences between reactant and product units

        correspo.reserve(products.size());

        {
            size_type num_units(reactants.size());

            size_type idx1 = 0;
            for (const_iterator i(products.begin()); i != products.end(); ++i, ++idx1)
            {
                size_type idx2 = 0;
                for (const_iterator j(reactants.begin()); j != reactants.end(); ++j, ++idx2)
                {
                    if (is_correspondent(*i, *j))
                    {
                        if (correspo.size() > idx1)
                        {
                            ;  //WARN: multiple correspondence found
                            assert(false);  // never get here
                        }
                        else if (std::find(correspo.begin(), correspo.end(), idx2)
                            != correspo.end())
                        {
                            ;  //WARN: multiple correspondence skipped
                            ;  // The reactant matches is already assigned to the other product
                        }
                        else
                        {
                            correspo.push_back(idx2);
                            break;
                        }
                    }
                }

                // If no reactant is found, create new one
                if (correspo.size() == idx1)
                {
                    correspo.push_back(num_units);
                    ++num_units;
                }
            }
        }

        // List reactants removed after reacting
        for (size_type i(0); i < reactants.size(); ++i)
        {
            if (std::find(correspo.begin(), correspo.end(), i) == correspo.end())
            {
                removed.push_back(i);
            }
        }

        return res;
    }

    typedef struct
    {
        std::vector<UnitSpecies> units;
        std::vector<unsigned int> groups;
        unsigned int num_groups;
    } unit_group_type;

    unit_group_type genunits(const operation_type& operations)
    {
        typedef std::vector<UnitSpecies>::size_type size_type;
        // typedef std::vector<UnitSpecies>::const_iterator const_iterator;

        if (itr_ != matchers_.end())
        {
            // Failed to match
            return unit_group_type();;
        }

        const std::vector<UnitSpecies>& products = operations.products;
        const std::vector<size_type>& correspo = operations.correspo;
        const std::vector<size_type>& removed = operations.removed;
        const size_type& reserved = operations.reserved;

        unit_group_type res = {};
        std::vector<UnitSpecies>& units = res.units;
        std::vector<unsigned int>& groups = res.groups;
        unsigned int& num_groups = res.num_groups;

        // 3. Concatenate units given as reactants

        const context_type ctx(context());

        int bond_stride = 0;

        for (typename element_container_type::const_iterator
            i(target_.begin()); i != target_.end(); ++i)
        {
            bond_stride += concatenate_units(units, (*i).units(), bond_stride);
        }

        // 4. Modify units

        utils::get_mapper_mf<std::string, std::string>::type bond_cache;

        {
            std::vector<std::pair<size_type, size_type> > priorities;
            {
                size_type next_idx = units.size();
                size_type idx = 0;
                for (std::vector<size_type>::const_iterator i(correspo.begin());
                    i != correspo.end(); ++i, ++idx)
                {
                    if ((*i) < reserved)
                    {
                        assert(ctx.iterators.size() > (*i));
                        priorities.push_back(std::make_pair(ctx.iterators[(*i)], idx));
                    }
                    else
                    {
                        priorities.push_back(std::make_pair(next_idx, idx));
                        ++next_idx;
                    }
                }
                std::sort(priorities.begin(), priorities.end());
            }

            for (std::vector<std::pair<size_type, size_type> >::const_iterator itr1(priorities.begin());
                itr1 != priorities.end(); ++itr1)
            {
                const UnitSpecies& op = products[(*itr1).second];
                const size_type& tgt = (*itr1).first;

                if (tgt >= units.size())
                {
                    // 4-1. Create a new unit
                    assert(tgt == units.size());
                    // tgt = units.size();
                    units.push_back(op);
                    if (rbex::is_named_wildcard(op.name()))
                    {
                        context_type::variable_container_type::const_iterator
                            itr(ctx.globals.find(op.name()));
                        if (itr == ctx.globals.end())
                        {
                            std::stringstream message;
                            message << "A named wildcard [" << op.name() << "] cannot be resolved.";
                            throw IllegalState(message.str());
                        }
                        else
                        {
                            units.back().set_name((*itr).second);
                        }
                    }
                }
                // else
                // {
                //     tgt = ctx.iterators[idx2];
                // }

                for (UnitSpecies::container_type::const_iterator i(op.begin());
                    i != op.end(); ++i)
                {
                    UnitSpecies::container_type::value_type&
                        site(units[tgt].at((*i).first));

                    // 4-2. Update sites
                    if ((*i).second.first == "")
                    {
                        ; // do nothing
                    }
                    else if (rbex::is_wildcard((*i).second.first))
                    {
                        if ((*i).second.first.size() != 1)
                        {
                            context_type::variable_container_type::const_iterator
                                itr(ctx.globals.find((*i).second.first));
                            if (itr == ctx.globals.end())
                            {
                                std::stringstream message;
                                message << "An invalid global name [" << (*i).second.first << "] was given.";
                                throw IllegalState(message.str());
                            }
                            else
                            {
                                site.second.first = (*itr).second;
                            }
                        }
                    }
                    else
                    {
                        site.second.first = (*i).second.first;
                    }

                    // 4-3. Update bonds
                    if ((*i).second.second == "")
                    {
                        // Just cut the bond
                        site.second.second = "";
                    }
                    else if (rbex::is_wildcard((*i).second.second))
                    {
                        if (rbex::is_named_wildcard((*i).second.second))
                        {
                            std::stringstream message;
                            message << "A named wildcard [" << (*i).second.second << "cannot be a bond."
                                << " Only a unnamed wildcard is accepted here.";
                            throw IllegalState(message.str());
                        }

                        ; // do nothing
                    }
                    else
                    {
                        // Make a new bond
                        utils::get_mapper_mf<std::string, std::string>::type::const_iterator
                            itr(bond_cache.find((*i).second.second));
                        if (itr != bond_cache.end())
                        {
                            site.second.second = (*itr).second;
                        }
                        else
                        {
                            ++bond_stride;
                            site.second.second = itos(bond_stride);
                            bond_cache.insert(std::make_pair((*i).second.second, site.second.second));
                        }
                    }
                }
            }
        }

        // 5. Remove units deleted

        if (removed.size() > 0)
        {
            for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
                i(removed.begin()); i != removed.end(); ++i)
            {
                units[ctx.iterators[(*i)]] = UnitSpecies();
            }
        }

        // The following is originally in group_units()

        // 6. Check connections between bonds and units

        utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type bondinfo;
        const unsigned int done = units.size();

        std::vector<std::vector<std::vector<UnitSpecies>::size_type> > connections;
        connections.resize(units.size());

        {
            unsigned int idx = 0;
            for (std::vector<UnitSpecies>::const_iterator i(units.begin());
                i != units.end(); ++i, ++idx)
            {
                for (UnitSpecies::container_type::const_iterator j((*i).begin());
                    j != (*i).end(); ++j)
                {
                    const std::string bond((*j).second.second);
                    if (bond == "")
                    {
                        continue;
                    }
                    else if (rbex::is_wildcard(bond))
                    {
                        std::stringstream ss;
                        ss << "A bond in a product is a wildcard ["
                            << (*i).name() << "(" << (*j).first << "^" << bond << ")].";
                        throw IllegalState(ss.str());
                    }

                    utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::iterator
                        itr(bondinfo.find(bond));
                    if (itr == bondinfo.end())
                    {
                        bondinfo[bond] = std::make_pair((*j).first, idx);
                    }
                    else
                    {
                        const unsigned int another = (*itr).second.second;
                        if (another == done)
                        {
                            std::stringstream ss;
                            ss << "A bond in a product is multiply connected ["
                                << (*i).name() << "(" << (*j).first << "^" << bond << ")].";
                            throw IllegalState(ss.str());
                        }

                        connections[idx].push_back(another);
                        connections[another].push_back(idx);
                        bondinfo[bond].second = done;  // This means the bond is already assigned.
                    }
                }
            }
        }

        if ((policy_ & ReactionRule::STRICT) && !(policy_ & (ReactionRule::DESTROY | ReactionRule::IMPLICIT)))
        {
            for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
                i(bondinfo.begin()); i != bondinfo.end(); ++i)
            {
                if ((*i).second.second != done)
                {
                    std::stringstream ss;
                    ss << "A bond in a product is not resolved [" << (*i).first << "].";
                    {
                        element_type sp;
                        for (std::vector<UnitSpecies>::const_iterator j(units.begin());
                            j != units.end(); ++j)
                        {
                            sp.add_unit((*j));
                        }
                        ss << " " << sp.serial();
                    }
                    throw IllegalState(ss.str());
                }
            }
        }

        // 7. Group units based on the connections

        const unsigned int notyet = units.size();
        groups.resize(units.size(), notyet);

        for (unsigned int idx(0); idx != units.size(); ++idx)
        {
            if (units[idx] == UnitSpecies())
            {
                continue;  // A removed unit
            }

            // std::cout << idx << " connects with ";
            // for (std::vector<std::vector<UnitSpecies>::size_type>::const_iterator
            //     i(connections[idx].begin()); i != connections[idx].end(); ++i)
            // {
            //     std::cout << *i;
            // }
            // std::cout << std::endl;

            num_groups = tag_units(
                groups, num_groups, idx, connections, notyet);
        }

        assert(num_groups <= notyet);

        // 8. Resolve open bonds based on the policy

        if (policy_ & ReactionRule::IMPLICIT)
        {
            for (std::vector<UnitSpecies>::iterator i(units.begin());
                i != units.end(); ++i)
            {
                for (UnitSpecies::container_type::size_type j(0);
                    j != (*i).num_sites(); ++j)
                {
                    UnitSpecies::container_type::value_type& site((*i).at(j));
                    const std::string bond(site.second.second);
                    if (bond != "" && !rbex::is_wildcard(bond) && bondinfo[bond].second != done)
                    {
                        site.second.second = "";
                    }
                }
            }
        }

        if (policy_ & ReactionRule::DESTROY)
        {
            for (utils::get_mapper_mf<std::string, std::pair<std::string, unsigned int> >::type::const_iterator
                i(bondinfo.begin()); i != bondinfo.end(); ++i)
            {
                if ((*i).second.second != done && units[(*i).second.second] != UnitSpecies())
                {
                    const unsigned int group_id = groups[(*i).second.second];

                    // assert(num_groups > 0);
                    // --num_groups;
                    std::vector<unsigned int>::iterator j1(groups.begin());
                    for (std::vector<UnitSpecies>::iterator j2(units.begin());
                        j2 != units.end(); ++j1, ++j2)
                    {
                        if ((*j1) == group_id)
                        {
                            (*j2) = UnitSpecies();
                            // (*j1) = notyet;
                        }
                    }
                }
            }
        }

        return res;
    }

    std::vector<element_type> generate()
    {
        const unit_group_type res = this->genunits(this->compile());
        return group_units<element_type>(res.units, res.groups, res.num_groups);
    }

    template <typename Generator>
    std::vector<typename Generator::return_type> gen(const element_container_type& reactants, const Generator& generator)
    {
        typedef std::vector<typename Generator::return_type> return_type;

        if (!this->match(reactants))
        {
            return return_type(0);
        }
        else if (pttrn_reactants_.size() == 0)
        {
            return return_type(
                1, generator(reactants, pttrn_products_));  // Zeroth-order reactions
            // return return_type(
            //     1, ReactionRule(Binder()(reactants), pttrn_products_, pttrn_.k()));  // Zeroth-order reactions
        }

        std::vector<std::vector<UnitSpecies> > candidates;
        const operation_type op = this->compile();

        return_type res;

        do
        {
            const unit_group_type _res = this->genunits(op);

            std::vector<std::vector<UnitSpecies> >::iterator i(std::find(candidates.begin(), candidates.end(), _res.units));
            if (i != candidates.end())
            {
                ; // (*i).set_k((*i).k() + rr.k());
            }
            else
            {
                candidates.push_back(_res.units);
                res.push_back(
                    generator(reactants, group_units<element_type>(_res.units, _res.groups, _res.num_groups)));
                // res.push_back(ReactionRule(Binder()(reactants), Binder()(group_units<element_type>(_res.units, _res.groups, _res.num_groups)), pttrn_.k()));
            }
        }
        while (this->next());

        return res;
    }

    const element_container_type& reactants() const
    {
        return target_;
    }

protected:

    const element_container_type& pttrn_reactants_;
    const element_container_type& pttrn_products_;
    const ReactionRule::policy_type& policy_;

    element_container_type target_;
    // std::vector<typename element_container_type::size_type> permutation_;
    std::vector<SpeciesExpressionMatcher> matchers_;
    std::vector<SpeciesExpressionMatcher>::iterator itr_;
};

struct Dummy
{
    std::vector<UnitSpecies> data;
    std::string _serial;

    explicit Dummy(Species const& sp)
        : data(sp.units()), _serial(sp.serial())
    {
        ; // do nothing
    }

    explicit Dummy(std::vector<UnitSpecies> const& units)
        : data(0), _serial("")
    {
        data.reserve(units.size());

        for (std::vector<UnitSpecies>::const_iterator i(units.begin());
            i != units.end(); ++i)
        {
            add_unit(*i);
        }
    }

    Dummy()
        : data()
    {
        ; // do nothing
    }

    const Species get() const
    {
        Species res;
        for (std::vector<UnitSpecies>::const_iterator i(data.begin());
            i != data.end(); ++i)
        {
            res.add_unit(*i);
        }
        return res;
    }

    const std::string serial() const
    {
        return _serial;
    }

    void add_unit(const UnitSpecies& usp)
    {
        data.push_back(usp);

        if (usp.name() == "")
        {
            throw NotSupported("UnitSpecies must have a name.");
        }
        else if (_serial != "")
        {
            _serial += "." + usp.serial();
        }
        else
        {
            _serial = usp.serial();
        }
    }

    const std::vector<UnitSpecies>& units() const
    {
        return data;
    }

    bool operator==(const Dummy& rhs) const
    {
        return serial() == rhs.serial();
    }

    bool operator!=(const Dummy& rhs) const
    {
        return serial() != rhs.serial();
    }

    bool operator<(const Dummy& rhs) const
    {
        return serial() < rhs.serial();
    }

    bool operator>(const Dummy& rhs) const
    {
        return serial() > rhs.serial();
    }
};

struct DummyReactionRule
{
    typedef DummyReactionRule return_type;

    ReactionRule data;
    std::vector<Dummy> _reactants, _products;

    DummyReactionRule()
        : data(), _reactants(), _products()
    {
        ;
    }

    DummyReactionRule(const ReactionRule& data)
        : data(data), _reactants(DummyBinder<Species, Dummy>::encode(data.reactants())), _products(DummyBinder<Species, Dummy>::encode(data.products()))
    {
        ;
    }

    DummyReactionRule(const ReactionRule& data, const std::vector<Dummy>& reactants, const std::vector<Dummy>& products)
        : data(data), _reactants(reactants), _products(products)
    {
        ;
    }

    const std::vector<Dummy>& reactants() const
    {
        return _reactants;
    }

    const std::vector<Dummy>& products() const
    {
        return _products;
    }

    ReactionRule get() const
    {
        ReactionRule res(DummyBinder<Species, Dummy>::decode(_reactants), DummyBinder<Species, Dummy>::decode(_products), data.k());
        res.set_policy(data.policy());
        return res;
    }

    std::vector<DummyReactionRule> generate(const std::vector<Dummy>& reactants) const
    {
        return ReactionRuleExpressionMatcher<Dummy>(_reactants, _products, data.policy()).gen(reactants, *this);
    }

    return_type operator()(const std::vector<Dummy>& reactants, const std::vector<Dummy>& products) const
    {
        return DummyReactionRule(data, reactants, products);
    }
};

DummyReactionRule format_reaction_rule(const DummyReactionRule& rr);

namespace rbex_extras
{

std::pair<std::vector<ReactionRule>, bool> apply_reaction_rules(
    std::vector<Species>& seeds,
    const std::vector<ReactionRule>& rules,
    const Integer max_itr,
    const std::map<Species, Integer>& max_stoich);

} // rbex_extras

} // ecell4

#endif /* ECELL4_CONTEXT_HPP */
