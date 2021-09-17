# -*- coding: utf-8 -*-
r"""
This file implements classes to deal with making changes to SimilaritySurfaces
and produce maps corresponding to those changes.
"""
######################################################################
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021      Pat Hooper
#
#  sage-flatsurf is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  sage-flatsurf is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sage-flatsurf. If not, see <https://www.gnu.org/licenses/>.
######################################################################

from .mappings import SurfaceMapping
from .polygon import ConvexPolygons, PolygonPosition, wedge_product, dot_product

from sage.groups.affine_gps.affine_group import AffineGroup
from sage.misc.cachefunc import cached_method
from sage.rings.ring import Ring
from sage.structure.sage_object import SageObject


class Partition(SageObject):
    r'''
    This class is represents a labeled disjoint collection of strictly convex
    polygons defined over a common base ring.

    Polygons should meet edge-to-edge.

    The `container()` method allows for determining which of the polygons
    contains a given tangent vector in the plane.

    A `Partition` can be constructed from a `base_ring` and a dictionary
    mapping labels to polygons.

    EXAMPLES::

        sage: from flatsurf.geometry.polygon import *
        sage: CP = ConvexPolygons(QQ)
        sage: poly3 = CP(vertices=[(2,2),(1,3),(0,2),(0,0)])
        sage: poly4 = CP(vertices=[(2,0),(2,2),(0,0)])
        sage: from flatsurf.geometry.mutation import Partition
        sage: part2 = Partition(QQ, {3:poly3, 4:poly4})
        sage: TestSuite(part2).run()
        sage: part2.container(vector(QQ,(1,2)), vector(QQ,(0,1)))
        3
    '''
    def __init__(self, base_ring, data, base_label = None):
        if not isinstance(base_ring, Ring):
            raise ValueError('base_ring should be a Ring')
        self._base_ring = base_ring
        self._data = data
        self._base_label = base_label

        internal_edges = {}
        external_edges = {}
        for label,p in self._data.items():
            ne = p.num_edges()
            for e in range(ne):
                vertex_pair = frozenset([p.vertex(e), p.vertex((e+1)%ne)])
                try:
                    other = external_edges.pop(vertex_pair)
                    internal_edges[vertex_pair] = (other, (label, e))
                except KeyError:
                    external_edges[vertex_pair] = (label, e)
        self._internal_edges = internal_edges
        self._external_edges = external_edges

    def base_label(self):
        return self._base_label

    def internal_edges(self, copy=True):
        r'''
        Return a dictionary mapping pairs of vertices (represented
        as a `frozenset` containing the two vertices), representing
        a line segment in the plane that appears as two edges of
        polygons in the partition to the pair of combinatorial
        edges `((label1, e1), (label2, e2))`.

        If `copy=True` (the default), a new dictionary is populated
        with the data. If `copy=False`, then the internal dictionary
        is returned. In this case, the user should be careful
        not to change the dictionary.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import *
            sage: CP = ConvexPolygons(QQ)
            sage: poly3 = CP(vertices=[(2,2),(1,3),(0,2),(0,0)])
            sage: poly4 = CP(vertices=[(2,0),(2,2),(0,0)])
            sage: from flatsurf.geometry.mutation import Partition
            sage: part2 = Partition(QQ, {3:poly3, 4:poly4})
            sage: part2.internal_edges()
            {frozenset({(0, 0), (2, 2)}): ((3, 3), (4, 1))}
        '''
        if copy:
            return dict(self._internal_edges)
        else:
            return self._internal_edges

    def external_edges(self, copy=True):
        r'''
        Return a dictionary mapping pairs of vertices (represented
        as a `frozenset` containing the two vertices), representing
        a line segment in the plane that appears as one edge of
        a polygon in the partition to the corresponding combinatorial
        edge pair `(label, edge)`.

        If `copy=True` (the default), a new dictionary is populated
        with the data. If `copy=False`, then the internal dictionary
        is returned. In this case, the user should be careful
        not to change the dictionary.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import *
            sage: CP = ConvexPolygons(QQ)
            sage: poly3 = CP(vertices=[(2,2),(1,3),(0,2),(0,0)])
            sage: poly4 = CP(vertices=[(2,0),(2,2),(0,0)])
            sage: from flatsurf.geometry.mutation import Partition
            sage: part2 = Partition(QQ, {3:poly3, 4:poly4})
            sage: d = part2.external_edges()
            sage: d[frozenset([vector(QQ,[0, 0],immutable=True),
            ....:              vector(QQ,[2, 0],immutable=True)])]
            (4, 2)
        '''
        if copy:
            return dict(self._external_edges)
        else:
            return self._external_edges

    def _test_polygons(self, **options):
        r"""
        Check that the polygons lie in `ConvexPolygons(base_ring)`.
        """
        tester = self._tester(**options)
        CP = ConvexPolygons(self.base_ring())
        for label, poly in self._data.items():
            tester.assertEqual(poly.parent(),CP,
                f'Polygon with label {label} doesn\'t have parent `ConvexPolygons(base_ring)`.')

    def base_ring(self):
        r'''
        Return the base ring.
        '''
        return self._base_ring

    def polygon(self, label):
        r'''
        Return the polygon associated to the given label.
        '''
        return self._data[label]

    def plot(self):
        r'''
        Plot the partition.
        '''
        it = iter(self._data.items())
        label,poly = next(it)
        plt = poly.plot()
        for label,poly in it:
            plt += poly.plot()
        return plt

    def labels(self):
        r'''
        Return the set of labels of polygons.
        '''
        return self._data.keys()

    def container(self, pt, vect):
        r'''
        Return the label of the polygon containing the tangent vector
        based at `pt` and with direction provided by `vect`.

        To be in a polygon, the tangent vector must be based in the interior
        of the polygon, or based on the boundary and pointed into the polygon,
        or tangent to the polygon and pointed in the counterclockwise
        direction around the polygon.

        EXAMPLES::

            sage: from flatsurf.geometry.polygon import *
            sage: CP = ConvexPolygons(QQ)
            sage: poly0 = CP(vertices=[(0,0),(2,0),(2,2),(1,2)])
            sage: poly1 = CP(vertices=[(1,2),(0,2),(0,0)])
            sage: poly2 = CP(vertices=[(0,2),(2,2),(1,3)])
            sage: from flatsurf.geometry.mutation import Partition
            sage: part1 = Partition(QQ, {0:poly0, 1:poly1, 2:poly2})
            sage: TestSuite(part1).run()
            sage: part1.container(vector(QQ,(1,1)),vector(QQ,(1,0)))
            0
            sage: part1.container(vector(QQ,(3/2,2)),vector(QQ,(1,0)))
            2
            sage: part1.container(vector(QQ,(3/2,2)),vector(QQ,(-1,0)))
            0
            sage: part1.container(vector(QQ,(1,2)),vector(QQ,(1,1)))
            2
            sage: part1.container(vector(QQ,(1,2)),vector(QQ,(1,0)))
            2
            sage: part1.container(vector(QQ,(1,2)),vector(QQ,(-1,0)))
            1
            sage: part1.container(vector(QQ,(1,2)),vector(QQ,(-1,-2)))
            0
        '''
        for label,poly in self._data.items():
            pos = poly.get_point_position(pt)
            if pos.is_inside():
                if pos.is_in_interior():
                    return label
                elif pos.is_in_edge_interior():
                    edge_vector = poly.edge(pos.get_edge())
                    wp = wedge_product(edge_vector,vect)
                    if wp > 0:
                        return label
                    if wp == 0:
                        if dot_product(edge_vector,vect) > 0:
                            return label
                elif pos.is_vertex():
                    v = pos.get_vertex()
                    edge_vector = poly.edge(v)
                    wp = wedge_product(edge_vector,vect)
                    if wp < 0:
                        continue
                    if wp == 0:
                        if dot_product(edge_vector,vect) < 0:
                            continue
                    edge_vector = poly.edge((v+poly.num_edges()-1)%poly.num_edges())
                    wp = wedge_product(edge_vector,vect)
                    if wp > 0:
                        return label
                else:
                    raise ValueError('Inside polygon but not in interior, in interior of an edge, or at a vertex. This is a bug.')
        # Searching for the container failed. Print some error:
        if vect == vect.parent().zero():
            raise ValueError('This method requires a non-zero `vect`.')
        raise ValueError(f'The pair ({pt}, {vect}) does not represent a tangent vector in the partition.')

    def __eq__(self, other):
        r'''
        Check that two partitions are identical.
        '''
        if not isinstance(other,Partition):
            raise ValueError('Can only compare to other partitions.')
        if not self.base_ring() == other.base_ring():
            return False
        if not self.base_label() == other.base_label():
            return False
        return self._data == other._data

class PartialMapToPlane(SageObject):
    r'''
    Members of this class define a map from a finite union
    of polygons in a SimilaritySurface to the plane. The maps
    are affine on polygons.

    The map can be applied to tangent vectors using the `to_plane`
    method.

    In applications, we assume the polygons have disjoint interiors
    in this case, we have an inverse map from the image back
    to the surface. This inverse is provided in the `to_surface`
    method.

    The main application of this is the case when you have two
    SimilaritySurfaces that differ by a natural polygonal
    operation such as gluing adjancent polygons along edges.
    Such a map can be described by two `PartialMapToPlane`
    with the same image from each surface. To apply the map,
    you push a vector forward to the plane from the domain
    under the first `PartialMapToPlane`, and then pull it back
    to the codomain under the second. This is realized through
    the `MapThroughPlane` class (which is obtainable through the
    `alter_surface` method).

    To construct a PartialMapToPlane, you need to pass the
    constructor a SimilaritySurface and a dictionary mapping
    labels to elements of the `AffineGroup` over the `base_ring`
    of the surface. (Similarities and 2x2 matrices are also okay
    since these can be converted into `AffineGroup` elements.)

    EXAMPLES::

        sage: from flatsurf import *
        sage: s = similarity_surfaces.example()
        sage: from sage.groups.affine_gps.affine_group import AffineGroup
        sage: AG = AffineGroup(2, QQ)
        sage: g0 = AG(matrix([[2, 0],[1, 1]]), vector([5,3]))
        sage: g1 = AG(matrix([[0, 2],[-1, 1]]), vector([9,5]))
        sage: from flatsurf.geometry.mutation import PartialMapToPlane
        sage: pm = PartialMapToPlane(s, {0:g0, 1:g1})
        sage: TestSuite(pm).run()
    '''
    def __init__(self, s, data, allow_mutable=False, ring=None, base_label = None):
        from flatsurf.geometry.similarity_surface import SimilaritySurface
        if not isinstance(s, SimilaritySurface):
            raise ValueError('s should be a SimilaritySurface')
        if not allow_mutable and s.is_mutable():
            raise ValueError('The surface must be immutable if `allow_mutable` is False.')
        self._s = s
        if ring:
            self._ring = ring
        else:
            self._ring = self._s.base_ring()
        self._setup(data, base_label)

    def _setup(self, data, base_label):
        r'''
        Initialize internal variables from data provided,
        where `data` is as in the `__init__` method.
        '''
        self._map = {}
        partition_data = {}
        AG = AffineGroup(2,self.base_ring())
        for label,ag in data.items():
            ag = AG(ag)
            self._map[label] = ag
            polygon = self._s.polygon(label).change_ring(self.base_ring())
            partition_data[label] = ag @ polygon
        self._part = Partition(self.base_ring(), partition_data, base_label=base_label)

    def surface(self):
        r'''
        Return the domain `SimilaritySurface` for this `PartialMapToPlane`.
        '''
        return self._s

    def base_ring(self):
        r'''
        Return the `base_ring` for the `SimilaritySurface` used.
        '''
        return self._ring

    def partition(self):
        r'''
        Return the `Partition` consisting of the images of the polygons
        in this `PartialMapToPlane`.
        '''
        return self._part

    def base_label(self):
        return self._part.base_label()

    def to_plane(self, tv):
        r'''
        Maps a tangent vector in the surface to a tangent vector in the plane.

        The tangent vector must be based at one of the polygons participating
        in this partial map. If this is not true then a `ValueError` is raised.

        The tangent vector is returned as a pair `(base_point, vector)`.
        '''
        label = tv.polygon_label()
        pt = tv.point()
        vect = tv.vector()
        try:
            ag = self._map[label]
        except KeyError:
            raise ValueError(f'The label `{label}` is not participating in this PartialMapToPlane.')
        return (ag(pt), ag.matrix()[:2,:2]*vect)

    def to_surface(self, pt, vect, ring = None):
        r'''
        Maps a tangent vector in the plane to a tangent vector on the surface.

        If the `ring` parameter is set, the returned vector will be in the tangent
        bundle of the surface defined over the given ring.

        The tangent vector in the plane is described by the a point in the plane,
        `pt`, and a vector in the plane `vect`.
        '''
        label = self._part.container(pt, vect)
        if ring == None or ring == self.base_ring():
            ag = ~self._map[label]
        else:
            ag = AffineGroup(2, ring)(~self._map[label])
        return self._s.tangent_vector(label, ag(pt), ag.matrix()[:2,:2]*vect, ring=ring)

    def __eq__(self, other):
        if not isinstance(other, PartialMapToPlane):
            raise ValueError('Can only compare with another PartialMapToPlane')
        if self._s != other._s:
            return False
        if self._map != other._map:
            return False
        return True

    def is_orientation_preserving(self):
        r'''
        Return `True` if all the affine maps used are orientation preserving,
        `False` if they are all orientation reversing, and `None` if both
        orientations occur.
        '''
        it = iter(self._map.items())
        label, ag = next(it)
        op = (ag.matrix().det() > 0)
        for label, ag in it:
            if not op == (ag.matrix().det() > 0):
                return None
        return op

    def is_affine_in_interior(self):
        r'''
        Return `True` if the map is continuous and locally affine in the interior
        of the domain. This method checks that edges that appear twice in the image
        (`intenal_edges()` of the partition) correspond to identified edges in the
        surface and that the map is affine along these edges.
        '''
        AG = AffineGroup(2, self.base_ring())
        for vertex_pair, ((l,e),(ll,ee)) in self._part.internal_edges().items():
            if self._s.opposite_edge(l,e) != (ll,ee):
                return False
            g = self._map[l]
            gg = self._map[ll]
            t1 = self._s.edge_transformation(l, e)
            if gg*AG(t1) != g:
                return False
        return True

    def _test_edge_vectors(self, **options):
        r"""
        We check that vectors pointed along the edges of the polygons
        in a counter-clockwise direction will be mapped to themselves
        when first pushed forward `to_plane` and then pulled back
        `to_surface`.
        """
        tester = self._tester(**options)
        s = self.surface()
        for label in self._map:
            p = s.polygon(label)
            for e in range(p.num_edges()):
                tv = s.tangent_vector(label, p.vertex(e), p.edge(e))
                pt, vect = self.to_plane(tv)
                tv2 = self.to_surface(pt, vect)
                tester.assertEqual(tv, tv2,
                    f'Failed edge vector test for edge {e} of polygon {label}.')

    def alter_surface(self, new_partition, in_place = False, partial_map = False, mapping = False):
        r'''
        Change the surface by replacing the polygons participating in this `PartialMapToPlane`
        with polygons from a new partition of the image.

        If `in_place` is `True`, then the changes are made to the surface from this
        `PartialMapToPlane`. In this case, this `PartialMapToPlane` must have been constructed
        with `allow_mutable=True`. This `PartialMapToPlane` is modified so that it maps polygons
        to the plane according to the new partition.

        The parameter `partial_map` may take one of three values: `True`, `False`, or the string
        `"mutable"`. If `True` or `"mutable"` a `PartialMapToPlane` will be returned. In this
        case the new surface can be recovered as `returned_value.surface()`. If `partial_map` is
        set to `True` and not `"mutable"`, this surface will not be mutable.

        If `mapping` is set to `True`, then a `MapThroughPlane` is returned whose domain is
        the surface associated to this `PartialMapToPlane` and whose codomain is the altered
        surface.

        EXAMPLES::

        In the following example, we take a surface built from two squares and cut one of the
        squares into three polygonal pieces, marking a new point at `(1/3, 2/3)`.

            sage: from flatsurf import *
            sage: S = SymmetricGroup(3)
            sage: right = S('(1,2)')
            sage: up = S('')
            sage: s = translation_surfaces.origami(right,up)
            sage: from sage.groups.affine_gps.affine_group import AffineGroup
            sage: AG = AffineGroup(2, QQ)
            sage: from flatsurf.geometry.mutation import PartialMapToPlane
            sage: pm = PartialMapToPlane(s, {1:AG.one()})
            sage: TestSuite(pm).run()
            sage: CP = ConvexPolygons(QQ)
            sage: poly0 = CP(vertices=[(0,0),(1,0),(1,1),(1/3,2/3)])
            sage: poly1 = CP(vertices=[(1/3,2/3),(1,1),(0,1)])
            sage: poly2 = CP(vertices=[(1/3,2/3),(0,1),(0,0)])
            sage: from flatsurf.geometry.mutation import Partition
            sage: new_partition = Partition(QQ, {0:poly0, 1:poly1, 2:poly2}, base_label=2)
            sage: s2 = pm.alter_surface(new_partition)
            sage: s2.num_polygons()
            4
            sage: TestSuite(s2).run()
            sage: mapping = pm.alter_surface(new_partition, mapping=True)
            sage: TestSuite(mapping).run()
            sage: mapping.is_cellular()
            False
            sage: mapping.is_invertible()
            True
            sage: mapping.is_locally_affine()
            True
            sage: mapping.is_orientation_preserving()
            True
            sage: s2.base_label()
            2
        '''
        if new_partition.external_edges().keys() != \
                self._part.external_edges().keys():
            raise ValueError('Partitions don\'t have the same external edges.')
        if new_partition.base_ring() != self.base_ring():
            raise ValueError('The partition must have the same base ring.')
        if in_place:
            if not self._s.is_mutable():
                raise ValueError('Can not alter immutable surface in place.')
            if mapping:
                raise ValueError('Can not alter surface in_place and also produce a mapping.')
            s2 = self._s
        else:
            s2 = self._s.copy(mutable=True)
        labels_to_be_removed = self.partition().labels()
        labels_to_be_added = new_partition.labels()
        edge_changes_forward = {}
        edge_changes_backward = {}
        for vertex_pair, (l,e) in self._part.external_edges().items():
            l_new, e_new = new_partition.external_edges()[vertex_pair]
            edge_changes_forward[(l,e)] = (True, l_new, e_new)
            # True indicates that this is a new edge.
        new_external_gluings = {}
        for vertex_pair, (l,e) in self._part.external_edges().items():
            ll,ee = self._s.opposite_edge(l,e)
            l_new, e_new = new_partition.external_edges()[vertex_pair]
            if ll in labels_to_be_removed:
                new_external_gluings[(l_new,e_new)] = edge_changes_forward[(ll,ee)]
            else:
                new_external_gluings[(l_new,e_new)] = (False, ll, ee)
                # False indicates that this is a new edge.
        us2 = s2.underlying_surface()
        if s2.base_label() in labels_to_be_removed:
            new_base_label = new_partition.base_label()
            if new_base_label == None:
                raise ValueError('When replacing the base polygon of a surface, you need to use a partition with a base_label.')
            us2.change_base_label(new_base_label)
        for label in labels_to_be_removed:
            if label in labels_to_be_added:
                us2.change_polygon(label, new_partition.polygon(label))
            else:
                us2.remove_polygon(label)
        altered_labels = {}
        for label in labels_to_be_added:
            if label not in labels_to_be_removed:
                try:
                    actual_label = us2.add_polygon(new_partition.polygon(label),
                                                   label=label)
                except ValueError:
                    actual_label = us2.add_polygon(new_partition.polygon(label))
                if actual_label != label:
                    altered_labels[label] = actual_label
        for (l,e),(new_edge,ll,ee) in new_external_gluings.items():
            if l in altered_labels:
                l = altered_labels[l]
            if new_edge and ll in altered_labels:
                ll = altered_labels[ll]
            us2.change_edge_gluing(l,e,ll,ee)
        for vertex_pair,((l,e),(ll,ee)) in new_partition.internal_edges().items():
            if l in altered_labels:
                l = altered_labels[l]
            if ll in altered_labels:
                ll = altered_labels[ll]
            us2.change_edge_gluing(l,e,ll,ee)
        # Deal with the base_label.

        # Construct data for the codomain PartialMapToPlane.
        AG = AffineGroup(2,self.base_ring())
        new_data = {}
        for label in labels_to_be_added:
            if label in altered_labels:
                label = altered_labels[label]
            new_data[label] = AG.one()
        # The data is now stored in new_data.
        if in_place:
            self._setup(new_data)
        if partial_map or mapping:
            if in_place:
                # partial_map must be true, but not mapping. This was checked above.
                return self
            else:
                AG = AffineGroup(2,self.base_ring())
                new_data = {}
                for label in labels_to_be_added:
                    if label in altered_labels:
                        label = altered_labels[label]
                    new_data[label] = AG.one()
                if partial_map:
                    if partial_map == 'mutable':
                        return PartialMapToPlane(s2, new_data, allow_mutable=True)
                    else:
                        us2.set_immutable()
                        return PartialMapToPlane(s2, new_data)
                else:
                    # return the mapping.
                    us2.set_immutable()
                    return MapThroughPlane(self, PartialMapToPlane(s2, new_data))
        else:
            # return the new surface
            return s2

class MapThroughPlane(SurfaceMapping):
    r"""
    This class makes it fairly easy to construct maps between surfaces via
    cut and paste operations using polygons.

    EXAMPLES::

    In the following example, we manually build an affine automorphism from the
    square torus to the hexagonal torus.

        sage: from flatsurf import *
        sage: s1 = translation_surfaces.square_torus()
        sage: s3 = translation_surfaces.veech_2n_gon(3)
        sage: ring = s3.base_ring()
        sage: sqrt3 = ring.gen()
        sage: sqrt3**2
        3

        sage: from flatsurf.geometry.mutation import PartialMapToPlane
        sage: from sage.groups.affine_gps.affine_group import AffineGroup
        sage: AG = AffineGroup(2, QQ)
        sage: partial_map = PartialMapToPlane(s1, {0:AG.one()})
        sage: CP = ConvexPolygons(QQ)
        sage: poly0 = CP(vertices=[(0,0),(1,0),(2/3,2/3),(0,1)])
        sage: poly1 = CP(vertices=[(1,0),(1,1),(2/3,2/3)])
        sage: poly2 = CP(vertices=[(0,1),(2/3,2/3),(1,1)])
        sage: from flatsurf.geometry.mutation import Partition
        sage: partition = Partition(QQ, {0:poly0, 1:poly1, 2:poly2}, base_label=1)

        sage: map1 = partial_map.alter_surface(partition, mapping=True)
        sage: TestSuite(map1).run()
        sage: s2 = map1.codomain()
        sage: s2.base_label()
        1
        sage: TestSuite(s2).run()
        sage: m = matrix(ring, [
        ....:    [     3/2,     0 ],
        ....:    [ sqrt3/2, sqrt3 ]
        ....: ])
        sage: AG2 = AffineGroup(2, ring)
        sage: from flatsurf.geometry.mutation import PartialMapToPlane
        sage: pm1 = PartialMapToPlane(s2, {
        ....:     0:AG2(m),
        ....:     1:AG2(m,vector( (-3/2, -sqrt3/2) )),
        ....:     2:AG2(m,vector( (0, -sqrt3) )),
        ....: }, ring = ring)
        sage: TestSuite(pm1).run()

        sage: pm2 = PartialMapToPlane(s3, {0: AG2.one()})
        sage: TestSuite(pm2).run()

        sage: from flatsurf.geometry.mutation import MapThroughPlane
        sage: map2 = MapThroughPlane(pm1, pm2)
        sage: TestSuite(map2).run()

        sage: phi = map2 * map1.change_ring(ring) # change required to make rings match
        sage: phi.domain()==s1 and phi.codomain()==s3
        True
        sage: TestSuite(phi).run()
        sage: phi.is_invertible()
        True
        sage: phi_inv = ~phi
        sage: TestSuite(phi_inv).run()
    """
    def __init__(self,
                 domain_partial_map_to_plane,
                 codomain_partial_map_to_plane):
        if not domain_partial_map_to_plane.base_ring() == codomain_partial_map_to_plane.base_ring():
            raise ValueError('Partial maps to the plane were defined over different base_rings.')
        self._map1 = domain_partial_map_to_plane
        self._map2 = codomain_partial_map_to_plane
        SurfaceMapping.__init__(self,
                                domain_partial_map_to_plane.surface(),
                                codomain_partial_map_to_plane.surface(),
                                ring = domain_partial_map_to_plane.base_ring()
                               )

    def support_contains(self, label):
        r'''
        Return if the associated polygon is in the support of the map.

        A polygon is in the support if the triple associated to a tangent
        vector in the polygon, `(label, pt, vect)`, is changed by the map.

        A return of `False` indicates that applying the map to
        `(label, pt, vect)` will lead to an equal return (with possibly
        the ring changing for `pt` and `vect`).

        A return of `True` suggests that a triple like this will be changed by
        the map, though we do not require this. By default this method returns
        `True`.
        '''
        return label in self._map1._map

    def domain_partial_map(self):
        r'''
        Return the domain `PartialMapToPlane` used to initialize this object.
        '''
        return self._map1

    def codomain_partial_map(self):
        r'''
        Return the codomain `PartialMapToPlane` used to initialize this object.
        '''
        return self._map2

    def push_vector_forward(self, tangent_vector, ring=None):
        r"""
        Applies the mapping to the provided vector. If the ring parameter
        is set, then the output will be a tangent vector in the provided ring.
        """
        try:
            pt, vect = self._map1.to_plane(tangent_vector)
        except ValueError:
            return self.codomain().tangent_vector(
                tangent_vector.polygon_label(),
                tangent_vector.point(),
                tangent_vector.vector(),
                ring=ring
            )
        return self._map2.to_surface(pt, vect, ring=ring)

    @cached_method
    def is_locally_affine(self):
        r"""
        Return `True` if the map is affine in local coordinates.

        Note that this requires that the map is affine at points in the
        interiors of edges.
        """
        if not self._map1.is_affine_in_interior():
            return False
        if not self._map2.is_affine_in_interior():
            return False

        # Now we need to check the boundary edges.
        AG = AffineGroup(2, self.base_ring())
        edge_to_transform = {}
        for vertex_pair, (l,e) in self._map1.partition().external_edges().items():
            image_l,image_e = self._map2.partition().external_edges()[vertex_pair]
            edge_to_transform[(l,e)] = (
                AG((~self._map2._map[image_l]) * self._map1._map[l]),
                image_l,
                image_e
            )
        for (l,e),(g,image_l,image_e) in edge_to_transform.items():
            ll,ee = self._map1._s.opposite_edge(l,e)
            if (ll,ee) in edge_to_transform:
                gg, image_ll, image_ee = edge_to_transform[(ll,ee)]
                if (image_ll, image_ee) != self._map2._s.opposite_edge(image_l,image_e):
                    # Not continuous
                    return False
            else:
                gg = AG.one()
                image_ll, image_ee = self._map2._s.opposite_edge(image_l,image_e)
            t1 = AG(self.domain().edge_transformation(l, e))
            t2 = AG(self.codomain().edge_transformation(image_l, image_e))
            if gg*t1 != t2*g:
                # Not affine along this edge
                return False
        return True

    @cached_method
    def is_cellular(self):
        r"""
        Return `True` if the map is cellular, i.e., polygons defining the
        surface in the domain are mapped bijectively to polygons in the
        codomain.

        We return `True` only if polygons in the two partitions are
        the same.
        """
        polys1 = {item[1] for item in self._map1._part._data.items()}
        polys2 = {item[1] for item in self._map2._part._data.items()}
        return polys1 == polys2

    @cached_method
    def is_orientation_preserving(self):
        r"""
        Return `True` if the mapping is orientation-preserving and
        `False` if not.

        We return True if both the partial maps to the plane are
        orientation preserving or both are orientation reversing.
        """
        op1 = self._map1.is_orientation_preserving()
        if op1 is None:
            return False
        op2 = self._map2.is_orientation_preserving()
        if op1 is None:
            return False
        return op1 == op2

    def is_invertible(self):
        r"""Return true if this mapping is invertible."""
        # Maps are not invertible by default
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.

        Raises a NotImplementedError by default.
        """
        return MapThroughPlane(self._map2, self._map1)
