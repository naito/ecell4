import collections
from cython cimport address
from cython.operator cimport dereference as deref, preincrement as inc
from ecell4.core cimport *


## ReactionInfo
cdef class ReactionInfo:
    """A class stores detailed information about a reaction in meso.

    ReactionInfo(t, reactants, products, coord)

    """

    def __init__(self, Real t, reactants, products, coord):
        """Constructor.

        Args:
          t (Real): A time when a reaction occurred.
          reactants (list): A list of reactants.
            Reactants are given as a ``Species``.
          products (list): A list of products.
            Products are given as a ``Species``.
          coord (int): A coordinate where a reaction occurred.

        """
        pass  #XXX: only used for doc string


    def __cinit__(self, Real t, reactants, products, Integer coord):
        cdef vector[Cpp_Species] reactants_
        cdef vector[Cpp_Species] products_

        for sp in reactants:
            reactants_.push_back(deref((<Species>sp).thisptr))
        for sp in products:
            products_.push_back(deref((<Species>sp).thisptr))

        self.thisptr = new Cpp_ReactionInfo(t, reactants_, products_, coord)

    def __dealloc__(self):
        del self.thisptr

    def t(self):
        """Return a time when a reaction occurred."""
        return self.thisptr.t()

    def coordinate(self):
        """Return a coordinate where a reaction occurred."""
        return self.thisptr.coordinate()

    def reactants(self):
        """Return a list of reactants

        Returns:
            list: A list of ``Species``.

        """
        cdef vector[Cpp_Species] species = self.thisptr.reactants()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def products(self):
        cdef vector[Cpp_Species] species = self.thisptr.products()
        """Return a list of products

        Returns:
            list: A list of ``Species``.

        """

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def __reduce__(self):
        return (ReactionInfo, (self.t(), self.reactants(), self.products(), self.coordinate()))

cdef ReactionInfo ReactionInfo_from_Cpp_ReactionInfo(Cpp_ReactionInfo* ri):
    cdef Cpp_ReactionInfo *new_obj = new Cpp_ReactionInfo(<Cpp_ReactionInfo> deref(ri))
    r = ReactionInfo(0, [], [], 0)
    del r.thisptr
    r.thisptr = new_obj
    return r

## MesoscopicWorld
#  a python wrapper for Cpp_MesoscopicWorld
cdef class MesoscopicWorld:
    """A class containing the properties of the meso world.

    MesoscopicWorld(edge_lengths=None, matrix_sizes=None, GSLRandomNumberGenerator rng=None)

    """

    def __init__(self, edge_lengths = None, matrix_sizes = None,
                 GSLRandomNumberGenerator rng = None):
        """Constructor.

        Parameters
        ----------
        edge_lengths : Real3, optional
            A size of the World.
        matrix_sizes : Integer3 or Real, optional
            A cell size of the space.
            The number of cells must be larger than 3, in principle.
            When given as a value, it is assumed to be a edge length of a cell.
        rng : GSLRandomNumberGenerator, optional
            A random number generator.

        """
        pass

    def __cinit__(self, edge_lengths = None,
        matrix_sizes = None, GSLRandomNumberGenerator rng = None):
        cdef string filename

        if edge_lengths is None:
            self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](new Cpp_MesoscopicWorld())
        elif matrix_sizes is None:
            if isinstance(edge_lengths, Real3):
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(deref((<Real3>edge_lengths).thisptr)))
            else:
                filename = tostring(edge_lengths)
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(filename))
        elif rng is None:
            if isinstance(matrix_sizes, Integer3):
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        deref((<Integer3>matrix_sizes).thisptr)))
            else:
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        <Real>matrix_sizes))
        else:
            if isinstance(matrix_sizes, Integer3):
                # XXX: GSLRandomNumberGenerator -> RandomNumberGenerator
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        deref((<Integer3>matrix_sizes).thisptr), deref(rng.thisptr)))
            else:
                self.thisptr = new shared_ptr[Cpp_MesoscopicWorld](
                    new Cpp_MesoscopicWorld(
                        deref((<Real3>edge_lengths).thisptr),
                        <Real>matrix_sizes, deref(rng.thisptr)))

    def __dealloc__(self):
        # XXX: Here, we release shared pointer,
        #      and if reference count to the MesoscopicWorld object,
        #      it will be released automatically.
        del self.thisptr

    def set_t(self, Real t):
        """set_t(t)

        Set the value of the time of the world.

        Parameters
        ----------
        t : Real
            The time of the world

        """
        self.thisptr.get().set_t(t)

    def t(self):
        """Return the time of the world."""
        return self.thisptr.get().t()

    def edge_lengths(self):
        """edge_lengths() -> Real3

        Return the edge lengths of the world.

        """
        cdef Cpp_Real3 lengths = self.thisptr.get().edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def actual_lengths(self):
        """actual_lengths() -> Real3

        Return the actual edge lengths of the world.
        Same as ``edge_lengths``.
        """
        cdef Cpp_Real3 lengths = self.thisptr.get().actual_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def matrix_sizes(self):
        """matrix_sizes() -> Integer3

        Return the number of subvolumes along axes.

        """
        cdef Cpp_Integer3 sizes = self.thisptr.get().matrix_sizes()
        return Integer3_from_Cpp_Integer3(address(sizes))

    def volume(self):
        """Return the volume of the world."""
        return self.thisptr.get().volume()

    def subvolume(self):
        """Return the subvolume of each cell."""
        return self.thisptr.get().subvolume()

    def get_value(self, Species sp):
        """get_value(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value(deref(sp.thisptr))

    def get_value_exact(self, Species sp):
        """get_value_exact(sp) -> Real

        Return the value (number) corresponding the given Species.

        Parameters
        ----------
        sp : Species
            a species whose value you require

        Returns
        -------
        Real:
            the value

        """
        return self.thisptr.get().get_value_exact(deref(sp.thisptr))

    def set_value(self, Species sp, Real value):
        """set_value(sp, value)

        Set the value of the given species.

        Parameters
        ----------
        sp : Species
            a species whose value you set
        value : Real
            a value set

        """
        self.thisptr.get().set_value(deref(sp.thisptr), value)

    def num_subvolumes(self, sp = None):
        """num_subvolumes(sp=None) -> Integer

        Return the number of subvolumes.

        Parameters
        ----------
        sp : Species, optional
            A species specifying a structure.
            When no species is given, return the total number of subvolumes.

        Returns
        -------
        Integer:
            The number of subvolumes.

        """
        if sp is None:
            return self.thisptr.get().num_subvolumes()
        else:
            return self.thisptr.get().num_subvolumes(deref((<Species>sp).thisptr))

    def subvolume_edge_lengths(self):
        """Return the edge lengths of a subvolume."""
        cdef Cpp_Real3 lengths = self.thisptr.get().subvolume_edge_lengths()
        return Real3_from_Cpp_Real3(address(lengths))

    def global2coord(self, Integer3 g):
        """global2coord(g) -> Integer

        Transform a global coordinate to a coordinate.

        Parameters
        ----------
        g : Integer3
            a global coordinate

        Returns
        -------
        Integer:
            A coordinate

        """
        return self.thisptr.get().global2coord(deref(g.thisptr))

    def coord2global(self, Integer c):
        """coord2global(coord) -> Integer3

        Transform a coordinate to a global coordinate.

        """
        cdef Cpp_Integer3 g = self.thisptr.get().coord2global(c)
        return Integer3_from_Cpp_Integer3(address(g))

    def position2global(self, Real3 pos):
        """position2global(pos) -> Integer3

        Transform a position to a global coordinate.

        Parameters
        ----------
        pos : Real3
            A position

        Returns
        -------
        Integer3:
            A global coordinate

        """
        cdef Cpp_Integer3 g = self.thisptr.get().position2global(deref(pos.thisptr))
        return Integer3_from_Cpp_Integer3(address(g))

    def position2coordinate(self, Real3 pos):
        """position2coordinate(pos) -> Integer

        Transform a position to a coordinate.

        Parameters
        ----------
        pos : Real3
            A position

        Returns
        -------
        Integer:
            A coordinate

        """
        return self.thisptr.get().position2coordinate(deref(pos.thisptr))

    def num_molecules(self, Species sp, c = None):
        """num_molecules(sp, c=None) -> Integer

        Return the number of molecules within the suggested subvolume.

        Parameters
        ----------
        sp : Species
            A species whose molecules you count
        c : Integer, optional
            A coordinate.

        Returns
        -------
        Integer:
            The number of molecules (of a given species)

        """
        if c is None:
            return self.thisptr.get().num_molecules(deref(sp.thisptr))
        elif isinstance(c, Integer3):
            return self.thisptr.get().num_molecules(deref(sp.thisptr), deref((<Integer3>c).thisptr))
        else:
            return self.thisptr.get().num_molecules(deref(sp.thisptr), <Integer>c)

    def num_molecules_exact(self, Species sp, c = None):
        """num_molecules_exact(sp, c=None) -> Integer

        Return the number of molecules within the suggested subvolume.

        Parameters
        ----------
        sp : Species
            The species of molecules to count
        c : Integer, optional
            A coordinate.

        Returns
        -------
        Integer:
            The number of molecules of a given species

        """
        if c is None:
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr))
        elif isinstance(c, Integer3):
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr), deref((<Integer3>c).thisptr))
        else:
            return self.thisptr.get().num_molecules_exact(deref(sp.thisptr), <Integer>c)

    def add_molecules(self, Species sp, Integer num, c = None):
        """add_molecules(sp, num, c=None)

        Add some molecules.

        Parameters
        ----------
        sp : Species
            a species of molecules to add
        num : Integer
            the number of molecules to add
        c : Integer or Shape, optional
            a coordinate or shape to add molecules on

        """
        if c is None:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num)
        elif isinstance(c, Integer3):
            self.thisptr.get().add_molecules(deref(sp.thisptr), num, deref((<Integer3>c).thisptr))
        elif hasattr(c, "as_base"):
            self.thisptr.get().add_molecules(
                deref(sp.thisptr), num, deref((<Shape>(c.as_base())).thisptr))
        else:
            self.thisptr.get().add_molecules(deref(sp.thisptr), num, <Integer>c)

    def remove_molecules(self, Species sp, Integer num, c = None):
        """remove_molecules(sp, num, c=None)

        Remove the molecules.

        Parameters
        ----------
        sp : Species
            a species whose molecules to remove
        num : Integer
            a number of molecules to be removed
        c : Integer, optional
            A coordinate.

        """
        if c is None:
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num)
        elif isinstance(c, Integer3):
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num, deref((<Integer3>c).thisptr))
        else:
            self.thisptr.get().remove_molecules(deref(sp.thisptr), num, <Integer>c)

    def add_structure(self, Species sp, shape):
        """add_structure(sp, shape)

        Add a structure.

        Parameters
        ----------
        sp : Species
            a species suggesting the shape.
        shape : Shape
            a shape of the structure.

        """
        self.thisptr.get().add_structure(
            deref(sp.thisptr), deref((<Shape>(shape.as_base())).thisptr))

    def get_volume(self, Species sp):
        """get_volume(sp) -> Real

        Return a volume of the given structure.

        Parameters
        ----------
        sp : Species
            A species for the target structure.

        Returns
        -------
        Real:
            A total volume of subvolumes belonging to the structure.

        """
        return self.thisptr.get().get_volume(deref(sp.thisptr))

    def has_structure(self, Species sp):
        """has_structure(sp) -> bool

        Return if the given structure is in the space or not.

        Parameters
        ----------
        sp : Species
            A species for the target structure.

        Returns
        -------
        bool:
            True if the given structure is in self.

        """
        return self.thisptr.get().has_structure(deref(sp.thisptr))

    def on_structure(self, Species sp, Integer3 g):
        """on_structure(sp, g) -> bool

        Check if the given species would be on the proper structure at the coordinate.

        Parameters
        ----------
        sp : Species
            a species scheduled to be placed
        g : Integer3
            a global coordinate pointing a subvolume

        Returns
        -------
        bool:
            True if it is on the proper structure.

        """
        return self.thisptr.get().on_structure(deref(sp.thisptr), deref(g.thisptr))

    def check_structure(self, Species sp, Integer3 g):
        """check_structure(sp, g) -> bool

        Check if the given subvolume is belonging to the structure.

        Parameters
        ----------
        sp : Species
            A species for the target structure.
        g : Integer3
            a global coordinate pointing a subvolume

        Returns
        -------
        bool:
            True if the subvolume is belonging to the structure.

        """
        return self.thisptr.get().check_structure(deref(sp.thisptr), deref(g.thisptr))

    def get_occupancy(self, Species sp, g):
        """get_occupancy(sp, g) -> Real

        Return the occupancy of the structure in the subvolume.

        Parameters
        ----------
        sp : Species
            A species for the target structure.
        g : Integer3
            a global coordinate pointing a subvolume

        Returns
        -------
        Real:
            The occupancy of the structure.
            As a default, return 1 if the subvolume overlaps with the structure,
            and return 0 otherwise.

        """
        if isinstance(g, Integer3):
            return self.thisptr.get().get_occupancy(
                deref(sp.thisptr), deref((<Integer3>g).thisptr))
        else:
            return self.thisptr.get().get_occupancy(deref(sp.thisptr), <Integer>g)

    def list_species(self):
        """Return a list of species."""
        cdef vector[Cpp_Species] species = self.thisptr.get().list_species()

        retval = []
        cdef vector[Cpp_Species].iterator it = species.begin()
        while it != species.end():
            retval.append(
                 Species_from_Cpp_Species(
                     <Cpp_Species*>(address(deref(it)))))
            inc(it)
        return retval

    def list_coordinates(self, Species sp):
        """list_coordinates(sp) -> [Integer]

        Return a list of coordinates of molecules belonging to the given species.

        Parameters
        ----------
        sp : Species
            A species of molecules.

        Returns
        -------
        list:
            A list of coordinates.

        """
        return self.thisptr.get().list_coordinates(deref(sp.thisptr))

    def list_coordinates_exact(self, Species sp):
        """list_coordinates_exact(sp) -> [Integer]

        Return a list of coordinates of molecules belonging to the given species.

        Parameters
        ----------
        sp : Species
            A species of molecules.

        Returns
        -------
        list:
            A list of coordinates.

        """
        return self.thisptr.get().list_coordinates_exact(deref(sp.thisptr))

    def new_particle(self, arg1, Real3 arg2=None):
        """new_particle(arg1, arg2=None) -> (ParticleID, Particle)

        Create a new particle.

        Parameters
        ----------
        arg1 : Particle
            A particle to be placed.

        or

        arg1 : Species
            A species of a particle
        arg2 : Real3
            A position to place a particle

        Returns
        -------
        tuple:
            A pair of ParticleID and Particle of a new particle

        """
        cdef pair[pair[Cpp_ParticleID, Cpp_Particle], bool] retval

        if arg2 is None:
            retval = self.thisptr.get().new_particle(deref((<Particle> arg1).thisptr))
        else:
            retval = self.thisptr.get().new_particle(deref((<Species> arg1).thisptr), deref(arg2.thisptr))
        return ((ParticleID_from_Cpp_ParticleID(address(retval.first.first)), Particle_from_Cpp_Particle(address(retval.first.second))), retval.second)

    def list_particles(self, Species sp = None):
        """list_particles(sp) -> [(ParticleID, Particle)]

        Return the list of particles.

        Parameters
        ----------
        sp : Species, optional
            The species of particles to list up
            If no species is given, return the whole list of particles.

        Returns
        -------
        list:
            The list of particles (of the given species)

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        if sp is None:
            particles = self.thisptr.get().list_particles()
        else:
            particles = self.thisptr.get().list_particles(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def list_particles_exact(self, Species sp):
        """list_particles_exact(sp) -> [(ParticleID, Particle)]

        Return the list of particles of a given species.

        Parameters
        ----------
        sp : Species
            The species of particles to list up

        Returns
        -------
        list:
            The list of particles of a given species

        """
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]] particles
        particles = self.thisptr.get().list_particles_exact(deref(sp.thisptr))

        retval = []
        cdef vector[pair[Cpp_ParticleID, Cpp_Particle]].iterator \
            it = particles.begin()
        while it != particles.end():
            retval.append(
                (ParticleID_from_Cpp_ParticleID(
                     <Cpp_ParticleID*>(address(deref(it).first))),
                 Particle_from_Cpp_Particle(
                     <Cpp_Particle*>(address(deref(it).second)))))
            inc(it)
        return retval

    def save(self, filename):
        """save(filename)

        Save the world to a file.

        Parameters
        ----------
        filename : str
            a filename to save to

        """
        self.thisptr.get().save(tostring(filename))

    def load(self, filename):
        """load(filename)

        Load the world from a file.

        Parameters
        ----------
        filename : str
            a filename to load from

        """
        self.thisptr.get().load(tostring(filename))

    def bind_to(self, m):
        """bind_to(m)

        Bind a model to the world

        Parameters
        ----------
        m : Model
            a model to bind

        """
        self.thisptr.get().bind_to(Cpp_Model_from_Model(m))

    def rng(self):
        """Return a random number generator object."""
        return GSLRandomNumberGenerator_from_Cpp_RandomNumberGenerator(
            self.thisptr.get().rng())

    def as_base(self):
        """Return self as a base class. Only for developmental use."""
        retval = Space()
        del retval.thisptr
        retval.thisptr = new shared_ptr[Cpp_Space](
            <shared_ptr[Cpp_Space]>deref(self.thisptr))
        return retval

cdef MesoscopicWorld MesoscopicWorld_from_Cpp_MesoscopicWorld(
    shared_ptr[Cpp_MesoscopicWorld] w):
    r = MesoscopicWorld(Real3(1, 1, 1), Integer3(1, 1, 1))
    r.thisptr.swap(w)
    return r

## MesoscopicSimulator
#  a python wrapper for Cpp_MesoscopicSimulator
cdef class MesoscopicSimulator:
    """ A class running the simulation with the meso algorithm.

    MesoscopicSimulator(m, w)

    """

    def __init__(self, m, MesoscopicWorld w=None):
        """MesoscopicSimulator(m, w)
        MesoscopicSimulator(w)

        Constructor.

        Parameters
        ----------
        m : Model
            A model
        w : MesoscopicWorld
            A world

        """
        pass

    def __cinit__(self, m, MesoscopicWorld w=None):
        if w is None:
            self.thisptr = new Cpp_MesoscopicSimulator(
                deref((<MesoscopicWorld>m).thisptr))
        else:
            self.thisptr = new Cpp_MesoscopicSimulator(
                Cpp_Model_from_Model(m), deref(w.thisptr))

    def __dealloc__(self):
        del self.thisptr

    def num_steps(self):
        """Return the number of steps."""
        return self.thisptr.num_steps()

    def step(self, upto = None):
        """step(upto=None) -> bool

        Step the simulation.

        Parameters
        ----------
        upto : Real, optional
            The time which to step the simulation up to

        Returns
        -------
        bool:
            True if the simulation did not reach the given time.
            When upto is not given, nothing will be returned.

        """
        if upto is None:
            self.thisptr.step()
        else:
            return self.thisptr.step(<Real> upto)

    def t(self):
        """Return the time."""
        return self.thisptr.t()

    def dt(self):
        """Return the step interval."""
        return self.thisptr.dt()

    def next_time(self):
        """Return the scheduled time for the next step."""
        return self.thisptr.next_time()

    def check_reaction(self):
        """Return if any reaction occurred at the last step, or not."""
        return self.thisptr.check_reaction()

    def last_reactions(self):
        """last_reactions() -> [(ReactionRule, ReactionInfo)]

        Return reactions occuring at the last step.

        Returns
        -------
        list:
            The list of reaction rules and infos.

        """
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]] reactions = self.thisptr.last_reactions()
        cdef vector[pair[Cpp_ReactionRule, Cpp_ReactionInfo]].iterator it = reactions.begin()
        retval = []
        while it != reactions.end():
            retval.append((
                ReactionRule_from_Cpp_ReactionRule(
                    <Cpp_ReactionRule*>(address(deref(it).first))),
                ReactionInfo_from_Cpp_ReactionInfo(
                    <Cpp_ReactionInfo*>(address(deref(it).second)))))
            inc(it)
        return retval

    def set_t(self, Real new_t):
        """set_t(t)

        Set the current time.

        Parameters
        ----------
        t : Real
            A current time.

        """
        self.thisptr.set_t(new_t)

    def set_dt(self, Real dt):
        """set_dt(dt)

        Set a step interval.

        Parameters
        ----------
        dt : Real
            A step interval

        """
        self.thisptr.set_dt(dt)

    def initialize(self):
        """Initialize the simulator."""
        self.thisptr.initialize()

    def model(self):
        """Return the model bound."""
        return Model_from_Cpp_Model(self.thisptr.model())

    def world(self):
        """Return the world bound."""
        return MesoscopicWorld_from_Cpp_MesoscopicWorld(self.thisptr.world())

    def run(self, Real duration, observers=None):
        """run(duration, observers)

        Run the simulation.

        Parameters
        ----------
        duration : Real
            A duration for running a simulation.
            A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            Observers

        """
        cdef vector[shared_ptr[Cpp_Observer]] tmp

        if observers is None:
            self.thisptr.run(duration)
        elif isinstance(observers, collections.Iterable):
            for obs in observers:
                tmp.push_back(deref((<Observer>(obs.as_base())).thisptr))
            self.thisptr.run(duration, tmp)
        else:
            self.thisptr.run(duration,
                deref((<Observer>(observers.as_base())).thisptr))

cdef MesoscopicSimulator MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
    Cpp_MesoscopicSimulator* s):
    r = MesoscopicSimulator(
        Model_from_Cpp_Model(s.model()),
        MesoscopicWorld_from_Cpp_MesoscopicWorld(s.world()))
    del r.thisptr
    r.thisptr = s
    return r

## MesoscopicFactory
#  a python wrapper for Cpp_MesoscopicFactory
cdef class MesoscopicFactory:
    """ A factory class creating a MesoscopicWorld instance and a MesoscopicSimulator instance.

    MesoscopicFactory(matrix_sizes=None, subvolume_length=None)

    """

    def __init__(self, matrix_sizes=None, subvolume_length=None):
        """Constructor.

        Parameters
        ----------
        matrix_sizes : Integer3, optional
            A cell size of the space.
            The number of cells must be larger than 3, in principle.
        subvolume_length : Real, optional
            The candidate for a edge length of a cell.
            If both matrix_sizes and subvolume_length are given,
            subvolume_length will be ignored.

        """
        pass

    def __cinit__(self, Integer3 matrix_sizes=None, subvolume_length=None):
        self.thisptr = new Cpp_MesoscopicFactory(
            Cpp_MesoscopicFactory.default_matrix_sizes() if matrix_sizes is None else deref(matrix_sizes.thisptr),
            Cpp_MesoscopicFactory.default_subvolume_length() if subvolume_length is None else <Real>subvolume_length)

    def __dealloc__(self):
        del self.thisptr

    def rng(self, GSLRandomNumberGenerator rng):
        """rng(GSLRandomNumberGenerator) -> MesoscopicFactory

        Set a random number generator, and return self.

        """
        cdef Cpp_MesoscopicFactory *ptr = self.thisptr.rng_ptr(deref(rng.thisptr))
        assert ptr == self.thisptr
        return self

    def create_world(self, arg1=None):
        """create_world(arg1=None) -> MesoscopicWorld

        Return a MesoscopicWorld instance.

        Parameters
        ----------
        arg1 : Real3
            The lengths of edges of a MesoscopicWorld created

        or

        arg1 : str
            The path of a HDF5 file for MesoscopicWorld

        Returns
        -------
        MesoscopicWorld:
            the created world

        """
        if arg1 is None:
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](
                    self.thisptr.create_world()))
        elif isinstance(arg1, Real3):
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](
                    self.thisptr.create_world(deref((<Real3>arg1).thisptr))))
        elif isinstance(arg1, str):
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](self.thisptr.create_world(<string>(arg1))))
        else:
            return MesoscopicWorld_from_Cpp_MesoscopicWorld(
                shared_ptr[Cpp_MesoscopicWorld](self.thisptr.create_world(
                    Cpp_Model_from_Model(arg1))))

    def create_simulator(self, arg1, MesoscopicWorld arg2=None):
        """create_simulator(arg1, arg2) -> MesoscopicSimulator

        Return a MesoscopicSimulator instance.

        Parameters
        ----------
        arg1 : MesoscopicWorld
            A world

        or

        arg1 : Model
            A simulation model
        arg2 : MesoscopicWorld
            A world

        Returns
        -------
        MesoscopicSimulator:
            the created simulator

        """
        if arg2 is None:
            return MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
                self.thisptr.create_simulator(deref((<MesoscopicWorld>arg1).thisptr)))
        else:
            return MesoscopicSimulator_from_Cpp_MesoscopicSimulator(
                self.thisptr.create_simulator(
                    Cpp_Model_from_Model(arg1), deref(arg2.thisptr)))
