#include "MolecularType.hpp"
#include "VacantType.hpp"
#include "LatticeSpace.hpp"

namespace ecell4
{

LatticeSpace::LatticeSpace(const Position3& edge_lengths,
        const Real& voxel_radius, const bool is_periodic) :
    voxel_radius_(voxel_radius), edge_lengths_(edge_lengths), t_(0),
    is_periodic_(is_periodic)
{
    vacant_ = new VacantType();
    border_ = new MolecularType(Species("Border", "0"));
    periodic_ = new MolecularType(Species("Periodic", "0"));

    set_lattice_properties(is_periodic_);
}

LatticeSpace::~LatticeSpace()
{
    delete vacant_;
    delete border_;
    delete periodic_;
}

/*
 * derived from SpatiocyteStepper::setLatticeProperties()
 */
void LatticeSpace::set_lattice_properties(const bool is_periodic)
{
    HCP_L = voxel_radius_/sqrt(3);
    HCP_X = voxel_radius_*sqrt(8.0/3); //Lx
    HCP_Y = voxel_radius_*sqrt(3); //Ly

    const Real lengthX = edge_lengths_[0];
    const Real lengthY = edge_lengths_[1];
    const Real lengthZ = edge_lengths_[2];

    row_size_ = (Integer)rint((lengthZ/2)/voxel_radius_) + 2;
    layer_size_ = (Integer)rint(lengthY/HCP_Y) + 2;
    col_size_ = (Integer)rint(lengthX/HCP_X) + 2;

    private_coordinate_type voxel_size(
            col_size_ * row_size_ * layer_size_);
    voxels_.reserve(voxel_size);
    for (private_coordinate_type coord(0); coord < voxel_size; ++coord)
    {
        Global global(private_coord2private_global(coord));
        if (global.col == 0 || global.col == col_size_-1 ||
                global.row == 0 || global.row == row_size_-1 ||
                global.layer == 0 || global.layer == layer_size_-1)
        {
            if (is_periodic)
                voxels_.push_back(periodic_);
            else
                voxels_.push_back(border_);
        }
        else
        {
            voxels_.push_back(vacant_);
        }
    }
}

Integer LatticeSpace::num_species() const
{
    return spmap_.size();
}

Integer LatticeSpace::num_molecules(const Species& sp) const
{
    return num_particles(sp) * 1;
}

const Position3& LatticeSpace::edge_lengths() const
{
    return edge_lengths_;
}

Integer LatticeSpace::num_particles() const
{
    Integer count(0);
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        count += (*itr).second.voxels().size();
    }

    return count;
}

Integer LatticeSpace::num_particles(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    Integer count(0);
    if (itr != spmap_.end())
    {
        count = (*itr).second.voxels().size();
    }
    return count;
}

bool LatticeSpace::has_particle(const ParticleID& pid) const
{
    return has_voxel(pid);
}

bool LatticeSpace::has_voxel(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return false;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return true;
            }
        }
    }

    return false;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles() const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }

        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::vector<std::pair<ParticleID, Particle> >
    LatticeSpace::list_particles(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Particle> > retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr != spmap_.end())
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            return retval;
        }

        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            retval.push_back(std::pair<ParticleID, Particle>(
                        (*vitr).second, particle_at((*vitr).first)));
        }
    }

    return retval;
}

std::pair<ParticleID, Voxel>
LatticeSpace::get_voxel(const ParticleID& pid) const
{
    for (spmap::const_iterator i(spmap_.begin()); i != spmap_.end(); ++i)
    {
        const MolecularType& mt((*i).second);
        typename MolecularType::container_type::const_iterator j(mt.find(pid));
        if (j != mt.end())
        {
            const coordinate_type coord(private2coord((*j).first));
            return std::make_pair(pid, Voxel((*i).first, coord, mt.radius(), mt.D()));
        }
    }

    throw NotFound("voxel not found.");
}

std::pair<ParticleID, Particle>
LatticeSpace::get_particle(const ParticleID& pid) const
{
    const Voxel v(get_voxel(pid).second);
    return std::make_pair(pid, Particle(
        v.species(), coordinate2position(v.coordinate()), v.radius(), v.D()));
}

bool LatticeSpace::update_particle(const ParticleID& pid, const Particle& p)
{
    return update_voxel_private(pid, Voxel(p.species(), position2private_coord(p.position()), p.radius(), p.D()));
}

/*
 * original methods
 */


std::vector<Species> LatticeSpace::list_species() const
{
    std::vector<Species> keys;
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        keys.push_back((*itr).first);
    }
    return keys;
}

std::vector<LatticeSpace::coordinate_type> LatticeSpace::list_coords(const Species& sp) const
{
    std::vector<coordinate_type> retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return retval;
    }

    const MolecularTypeBase* mt(&((*itr).second));

    for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
            itr != mt->end(); ++itr)
    {
        retval.push_back(private2coord((*itr).first));
    }

    return retval;
}

std::vector<std::pair<ParticleID, Voxel> >
LatticeSpace::list_voxels(const Species& sp) const
{
    std::vector<std::pair<ParticleID, Voxel> > retval;
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return retval;
    }

    const MolecularTypeBase* mt(&((*itr).second));
    for (MolecularTypeBase::container_type::const_iterator itr(mt->begin());
        itr != mt->end(); ++itr)
    {
        retval.push_back(std::make_pair(
            (*itr).second,
            Voxel(sp, private2coord((*itr).first), mt->radius(), mt->D())));
    }
    return retval;
}

std::pair<LatticeSpace::spmap::iterator, bool>
LatticeSpace::__get_molecular_type(const Voxel& v)
{
    LatticeSpace::spmap::iterator itr(spmap_.find(v.species()));
    if (itr != spmap_.end())
    {
        return std::make_pair(itr, false);
    }

    MolecularType mt(v.species(), v.radius(), v.D());
    std::pair<LatticeSpace::spmap::iterator, bool> retval(
        spmap_.insert(LatticeSpace::spmap::value_type(v.species(), mt)));
    if (!retval.second)
    {
        throw "insert error.";
    }
    return retval;
}

MolecularTypeBase* LatticeSpace::find_molecular_type(const Species& sp)
{
    LatticeSpace::spmap::iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        throw NotFound("MolecularType nof found.");
    }
    return &((*itr).second);
}

MolecularTypeBase* LatticeSpace::get_molecular_type(const Voxel& v)
{
    return &((*(__get_molecular_type(v).first)).second);
}

// bool LatticeSpace::register_species(const Species& sp)
// {
//     return __get_molecular_type(sp).second;
// }

/*
 * Coordinate transformations
 */

LatticeSpace::coordinate_type LatticeSpace::global2coord(const Global& global) const
{
    return global2coord_(global, col_size(), row_size(), layer_size());
}

LatticeSpace::private_coordinate_type LatticeSpace::global2private_coord(
        const Global& global) const
{
    Global modified_global;
    modified_global.col = global.col + 1;
    modified_global.row = global.row + 1;
    modified_global.layer = global.layer + 1;
    return global2coord_(modified_global, col_size_, row_size_, layer_size_);
}

const Global LatticeSpace::coord2global(coordinate_type coord) const
{
    return coord2global_(coord, col_size(), row_size(), layer_size());
}

const Global LatticeSpace::private_coord2global(
        private_coordinate_type private_coord) const
{
    Global retval(private_coord2private_global(private_coord));
    retval.col--;
    retval.row--;
    retval.layer--;
    return retval;
}

/*
 * Protected functions
 */

LatticeSpace::private_coordinate_type LatticeSpace::get_coord(const ParticleID& pid) const
{
    for (spmap::const_iterator itr(spmap_.begin());
            itr != spmap_.end(); ++itr)
    {
        const MolecularTypeBase& mt((*itr).second);
        if (mt.is_vacant())
        {
            continue;
        }
        for (MolecularTypeBase::container_type::const_iterator vitr(mt.begin());
                vitr != mt.end(); ++vitr)
        {
            if ((*vitr).second == pid)
            {
                return (*vitr).first;
            }
        }
    }
    return -1;
}

MolecularTypeBase* LatticeSpace::get_molecular_type(private_coordinate_type coord) const
{
    // return voxels_.at(coord2private(coord));
    return voxels_.at(coord);
}

bool LatticeSpace::has_species(const Species& sp) const
{
    return spmap_.find(sp) != spmap_.end();
}

bool LatticeSpace::remove_particle(const ParticleID& pid)
{
    return remove_voxel(pid);
}

bool LatticeSpace::remove_voxel(const ParticleID& pid)
{
    for (spmap::iterator i(spmap_.begin()); i != spmap_.end(); ++i)
    {
        MolecularType& mt((*i).second);
        typename MolecularType::container_type::const_iterator j(mt.find(pid));
        if (j != mt.end())
        {
            const private_coordinate_type coord((*j).first);
            if (mt.removeVoxel(coord))
            {
                voxel_container::iterator itr(voxels_.begin() + coord);
                (*itr) = vacant_;
                return true;
            }
            return false;
        }
    }
    return false;
}

bool LatticeSpace::remove_molecule(const private_coordinate_type coord)
{
    voxel_container::iterator itr(voxels_.begin() + coord);
    MolecularTypeBase* mt(*itr);
    if (mt->is_vacant())
    {
        return false;
    }
    if (mt->removeVoxel(coord))
    {
        (*itr) = vacant_;
        return true;
    }
    return false;
}

bool LatticeSpace::move(coordinate_type from, coordinate_type to)
{
    const coordinate_type private_from(coord2private(from));
    const coordinate_type private_to(coord2private(to));
    return move_(private_from, private_to).second;
}

std::pair<LatticeSpace::coordinate_type, bool> LatticeSpace::move_to_neighbor(
        private_coordinate_type coord, Integer nrand)
{
    const private_coordinate_type neighbor(get_neighbor(coord, nrand));
    return move_(coord, neighbor);
}

std::pair<LatticeSpace::private_coordinate_type, bool> LatticeSpace::move_to_neighbor(MolecularTypeBase::particle_info& info,
        Integer nrand)
{
    const coordinate_type neighbor(get_neighbor(info.first, nrand));
    return move_(info, neighbor);
}

std::pair<LatticeSpace::private_coordinate_type, bool> LatticeSpace::move_(
        private_coordinate_type private_from, private_coordinate_type private_to)
{
    if (private_from == private_to)
    {
        std::cerr << " from == to ";
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }

    MolecularTypeBase* from_mt(voxels_.at(private_from));
    if (from_mt->is_vacant())
    {
        return std::pair<private_coordinate_type, bool>(private_from, true);
    }

    MolecularTypeBase* to_mt(voxels_.at(private_to));

    if (to_mt == border_)
    {
        std::cerr << " to_mt is border ";
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }
    else if (to_mt == periodic_)
    {
        private_to = apply_boundary_(private_to);
        to_mt = voxels_.at(private_to);
    }

    if (!to_mt->is_vacant())
    {
        std::cerr << " to_mt is " << to_mt->species().serial() << " ";
        return std::pair<private_coordinate_type, bool>(private_to, false);
    }

    MolecularTypeBase::container_type::iterator itr(from_mt->find(private_from));
    (*itr).first = private_to;
    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = vacant_;
    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

std::pair<LatticeSpace::private_coordinate_type, bool> LatticeSpace::move_(
        MolecularTypeBase::particle_info& info, private_coordinate_type private_to)
{
    const private_coordinate_type private_from(info.first);
    if (private_from == private_to)
    {
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }

    MolecularTypeBase* from_mt(voxels_.at(private_from));
    if (from_mt->is_vacant())
    {
        return std::pair<private_coordinate_type, bool>(private_from, true);
    }

    MolecularTypeBase* to_mt(voxels_.at(private_to));

    if (to_mt == border_)
    {
        return std::pair<private_coordinate_type, bool>(private_from, false);
    }
    else if (to_mt == periodic_)
    {
        private_to = apply_boundary_(private_to);
        to_mt = voxels_.at(private_to);
    }

    if (!to_mt->is_vacant())
    {
        return std::pair<private_coordinate_type, bool>(private_to, false);
    }

    info.first = private_to;
    voxel_container::iterator from_itr(voxels_.begin() + private_from);
    (*from_itr) = vacant_;
    voxel_container::iterator to_itr(voxels_.begin() + private_to);
    (*to_itr) = from_mt;

    return std::pair<private_coordinate_type, bool>(private_to, true);
}

LatticeSpace::private_coordinate_type LatticeSpace::get_neighbor(
        private_coordinate_type private_coord, Integer nrand) const
{
    const Integer NUM_COLROW(col_size_ * row_size_);
    const Integer NUM_ROW(row_size_);
    const bool odd_col((private_coord % NUM_COLROW / NUM_ROW) & 1);
    const bool odd_lay((private_coord / NUM_COLROW) & 1);

    switch(nrand)
    {
    case 1:
        return private_coord+1;
    case 2:
        return private_coord+(odd_col^odd_lay)-NUM_ROW-1;
    case 3:
        return private_coord+(odd_col^odd_lay)-NUM_ROW;
    case 4:
        return private_coord+(odd_col^odd_lay)+NUM_ROW-1;
    case 5:
        return private_coord+(odd_col^odd_lay)+NUM_ROW;
    case 6:
        return private_coord+(2*odd_col-1)*NUM_COLROW-NUM_ROW;
    case 7:
        return private_coord-(2*odd_col-1)*NUM_COLROW+NUM_ROW;
    case 8:
        return private_coord+(odd_col^odd_lay)-NUM_ROW-1;
    case 9:
        return private_coord+(odd_col^odd_lay)-NUM_ROW;
    case 10:
        return private_coord+(odd_col^odd_lay)+NUM_ROW-1;
    case 11:
        return private_coord+(odd_col^odd_lay)+NUM_ROW;
    }
    return private_coord-1;
}

const Particle LatticeSpace::particle_at(private_coordinate_type coord) const
{
    const MolecularTypeBase* ptr_mt(voxels_.at(coord));
    // const Species& sp = ptr_mt->species();
    // const Position3& pos = coordinate2position(coord);
    // const Real& radius = 0;
    // const Real& D = 0;
    // Particle particle(sp, pos, radius, D);
    // return particle;
    return Particle(
        ptr_mt->species(), coordinate2position(coord), ptr_mt->radius(), ptr_mt->D());
}

bool LatticeSpace::is_in_range(coordinate_type coord) const
{
    return coord >= 0 && coord < row_size() * layer_size() * col_size();
}

bool LatticeSpace::is_in_range_private(private_coordinate_type coord) const
{
    return coord >= 0 && coord < row_size_ * col_size_ * layer_size_;
}

/*
 * Coordinate transformations
 */

LatticeSpace::coordinate_type LatticeSpace::global2coord_(const Global& global,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    return global.row +
        row_size * (global.col) +
        row_size * col_size * (global.layer);
}

const Global LatticeSpace::coord2global_(coordinate_type coord,
        Integer col_size, Integer row_size, Integer layer_size) const
{
    /*
    Global retval;
    retval.col = coord / (row_size * layer_size);
    retval.layer = (coord % (row_size * layer_size)) / row_size;
    retval.row = (coord % (row_size * layer_size)) % row_size;
    return retval;
    */
    Global retval;
    const Integer NUM_COLROW(row_size * col_size);
    const Integer LAYER(coord / NUM_COLROW);
    const Integer SURPLUS(coord - LAYER * NUM_COLROW);
    const Integer COL(SURPLUS / row_size);
    retval.col = COL;
    retval.layer = LAYER;
    retval.row = SURPLUS - COL * row_size;
    return retval;
}

const Global LatticeSpace::private_coord2private_global(
        const private_coordinate_type private_coord) const
{
    return coord2global_(private_coord, col_size_, row_size_, layer_size_);
}

const LatticeSpace::private_coordinate_type LatticeSpace::private_global2private_coord(
        const Global private_global) const
{
    return global2coord_(private_global, col_size_, row_size_, layer_size_);
}

const Position3 LatticeSpace::global2position(const Global& global) const
{
    //the center point of a voxel
    Position3 position;
    position[0] = global.col * HCP_X;
    position[1] = (global.col % 2) * HCP_L + HCP_Y * global.layer;
    position[2] = (global.row * 2 + (global.layer + global.col) % 2)
        * voxel_radius_;
    return position;
}

const Global LatticeSpace::position2global(const Position3& pos) const
{
    Global global;
    global.col = (Integer)(pos[0] / HCP_X);
    global.layer = (Integer)((pos[1] - (global.col % 2) * HCP_L) / HCP_Y);
    global.row = (Integer)((pos[2] - (((global.layer + global.col) % 2) *
                    voxel_radius_)) /
            voxel_radius_ / 2);
    return global;
}

const Position3 LatticeSpace::coordinate2position(coordinate_type coord) const
{
    return global2position(coord2global(coord));
}

LatticeSpace::coordinate_type LatticeSpace::position2coordinate(const Position3& pos) const
{
    return global2coord(position2global(pos));
}

LatticeSpace::private_coordinate_type LatticeSpace::position2private_coord(const Position3& pos) const
{
    return global2private_coord(position2global(pos));
}

LatticeSpace::private_coordinate_type LatticeSpace::coord2private(coordinate_type coord) const
{
    return global2private_coord(coord2global(coord));
}

LatticeSpace::coordinate_type LatticeSpace::private2coord(private_coordinate_type private_coord) const
{
    return global2coord(private_coord2global(private_coord));
}

LatticeSpace::private_coordinate_type LatticeSpace::apply_boundary_(const private_coordinate_type& private_coord) const
{
    Global global(private_coord2private_global(private_coord));

    global.col = (global.col - 1) % col_size();
    global.row = (global.row - 1) % row_size();
    global.layer = (global.layer - 1) % layer_size();

    global.col = global.col < 0 ? global.col + col_size() + 1 : global.col + 1;
    global.row = global.row < 0 ? global.row + row_size() + 1 : global.row + 1;
    global.layer = global.layer < 0 ? global.layer + layer_size() + 1 : global.layer + 1;

    return private_global2private_coord(global);
}

Integer LatticeSpace::num_voxels(const Species& sp) const
{
    spmap::const_iterator itr(spmap_.find(sp));
    if (itr == spmap_.end())
    {
        return 0;
    }
    const MolecularTypeBase* mt(&((*itr).second));
    return mt->size();
}

Integer LatticeSpace::num_voxels() const
{
    Integer retval(0);
    for (spmap::const_iterator itr(spmap_.begin()); itr != spmap_.end(); ++itr)
    {
        retval += (*itr).second.size();
    }
    return retval;
}

bool LatticeSpace::update_voxel(const ParticleID& pid, const Voxel& v)
{
    return update_voxel_private(pid,
        Voxel(v.species(), coord2private(v.coordinate()), v.radius(), v.D()));
}

bool LatticeSpace::update_voxel_private(const Voxel& v)
{
    const LatticeSpace::private_coordinate_type coord(v.coordinate());
    MolecularTypeBase* src_mt(voxels_.at(coord));
    if (src_mt->is_vacant())
    {
        return false;
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    if (new_mt->is_vacant())
    {
        return false; // Vacant has no ParticleID. Call remove_voxel.
    }

    MolecularTypeBase::particle_info info(*(src_mt->find(coord)));
    src_mt->removeVoxel(coord);
    new_mt->addVoxel(info);
    voxel_container::iterator itr(voxels_.begin() + coord);
    (*itr) = new_mt;
    return true;
}

bool LatticeSpace::update_voxel_private(const ParticleID& pid, const Voxel& v)
{
    const LatticeSpace::private_coordinate_type to_coord(v.coordinate());
    if (!is_in_range_private(to_coord))
    {
        return false;
    }

    MolecularTypeBase* dest_mt(get_molecular_type(to_coord));
    if (!dest_mt->is_vacant())
    {
        if (dest_mt->species() == periodic_->species()
            || dest_mt->species() == border_->species())
        {
            throw NotSupported("The coordinate points a boundary.");
        }

        MolecularTypeBase::container_type::const_iterator
            dest_itr(dest_mt->find(to_coord));
        if (dest_itr == dest_mt->end())
        {
            throw IllegalState(
                "MolecularTypaBase [" + dest_mt->species().serial()
                + "] doesn't contain a coordinate.");
        }
        else if ((*dest_itr).second != pid)
        {
            return false;
        }
    }

    MolecularTypeBase* new_mt(get_molecular_type(v)); //XXX: need MoleculeInfo
    if (new_mt->is_vacant())
    {
        return false; // Vacant has no ParticleID.
    }

    const LatticeSpace::private_coordinate_type from_coord(get_coord(pid));
    if (from_coord != -1)
    {
        MolecularTypeBase* src_mt(voxels_.at(from_coord));
        src_mt->removeVoxel(from_coord);
        voxel_container::iterator itr(voxels_.begin() + from_coord);
        voxels_.erase(itr);
        voxels_.insert(itr, vacant_);
    }

    new_mt->addVoxel(MolecularTypeBase::particle_info(to_coord, pid));
    voxel_container::iterator itr(voxels_.begin() + to_coord);
    (*itr) = new_mt;
    return true;
}

} // ecell4