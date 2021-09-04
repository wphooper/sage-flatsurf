# -*- coding: utf-8 -*-
#*********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2016-2021 Pat Hooper
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
#*********************************************************************
r"""Mappings between translation surfaces."""

from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip
from six import iteritems

from flatsurf.geometry.polygon import ConvexPolygons, wedge_product
from flatsurf.geometry.surface import Surface, Surface_list, Surface_dict, ExtraLabel
from flatsurf.geometry.similarity_surface import SimilaritySurface

from sage.groups.affine_gps.affine_group import AffineGroup
from sage.matrix.constructor import matrix
from sage.rings.infinity import Infinity
from sage.structure.sage_object import SageObject

class SurfaceMapping(SageObject):
    r"""Abstract class for any mapping between surfaces."""

    def __init__(self, domain, codomain):
        self._domain=domain
        self._codomain=codomain

    def domain(self):
        r"""
        Return the domain of the mapping.
        """
        return self._domain

    def codomain(self):
        r"""
        Return the range of the mapping.
        """
        return self._codomain

    def is_invertible(self):
        r"""Return true if this mapping is invertible."""
        # Maps are not invertible by default
        return False

    def is_locally_affine(self):
        r"""
        Return `True` if the map is affine in local coordinates.

        Note that this requires that the map is affine at points in the
        interiors of edges.
        """
        # Maps are not locally affine by default
        return False

    def is_cellular(self):
        r"""
        Return `True` if the map is cellular, i.e., polygons defining the
        surface in the domain are mapped bijectively to polygons in the
        codomain.
        """
        # Maps are not cellular by default.
        return False

    def is_orientation_preserving(self):
        r"""
        Return `True` if the mapping is orientation-preserving and
        `False` if not.

        By default we return `True`.
        """
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.

        Raises a NotImplementedError by default.
        """
        raise NotImplementedError('__invert__ has not been implemented.')

    def push_vector_forward(self,tangent_vector, ring=None):
        r"""
        Applies the mapping to the provided vector. If the ring parameter
        is set, then the output will be a tangent vector in the provided ring.
        """
        raise NotImplementedError

    def __mul__(self,other):
        # Compose SurfaceMappings
        return SurfaceMappingComposition(other,self)

    def __rmul__(self,other):
        return SurfaceMappingComposition(self,other)

    def __call__(self, x, ring=None):
        from flatsurf.geometry.tangent_bundle import SimilaritySurfaceTangentVector
        if isinstance(x, SimilaritySurfaceTangentVector):
            return self.push_vector_forward(x, ring=ring)

    def __eq__(self, other):
        if self is other:
            return True
        if self.domain() != other.domain():
            return False
        if self.codomain() != other.codomain():
            return False
        l = self._domain.base_label()
        p = self._domain.polygon(l)
        pt = p.vertex(0)
        v1 = p.edge(0)
        tv1_self = self._domain.tangent_vector(l, pt, v1)
        tv1_other = other._domain.tangent_vector(l, pt, v1)
        if self.push_vector_forward(tv1_self) != other.push_vector_forward(tv1_other):
            return False
        v2 = - p.edge(p.num_edges()-1)
        tv2_self = self._domain.tangent_vector(l, pt, v2)
        tv2_other = other._domain.tangent_vector(l, pt, v2)
        if self.push_vector_forward(tv2_self) != other.push_vector_forward(tv2_other):
            return False
        if matrix([v1,v2]).det() == 0:
            raise ValueError('Base polygon of self.domain() is not strictly convex!')
        if self.is_locally_affine() and other.is_locally_affine():
            # If two tangent vectors based at the same point are linear
            # independent and have the same image under locally affine
            # maps, then the maps are the same!
            return True
        raise NotImplemented('Comparison of mappings is not fully implemented unless both maps are locally affine.')

    def __hash__(self):
        if self.is_locally_affine():
            l = self._domain.base_label()
            p = self._domain.polygon(l)
            pt = p.vertex(0)
            v1 = p.edge(0)
            iv1 = self.push_vector_forward(
                self._domain.tangent_vector(l, pt, v1)
            )
            v2 = - p.edge(p.num_edges()-1)
            if matrix([v1,v2]).det() == 0:
                raise ValueError('Base polygon of domain() is not strictly convex!')
            iv2 = self.push_vector_forward(
                self._domain.tangent_vector(l, pt, v2)
            )
            return 23 * hash(self.domain()) + 47 * hash(self.codomain()) + \
                173 * hash(iv1) + 59 * hash(iv2)
        else:
            raise NotImplemented('Hash is only implemented for locally affine maps.')


class SurfaceMappingComposition(SurfaceMapping):
    r"""
    Composition of mappings.
    """

    def __init__(self, *mappings):
        r"""
        Represent the composition of mappings.

        Mappings must be listed in the order in which they will be applied.
        """
        if len(mappings) < 2:
            raise ValueError('Need at least two mappings to compose.')
        for i in range(1, len(mappings)):
            if mappings[i-1].codomain() != mappings[i].domain():
                raise ValueError(f'The codomain of map {i-1} does not match the domain of map {i}.')
        self._m = []
        for m in mappings:
            if isinstance(m, SurfaceMappingComposition):
                self._m += m.factors()
            else:
                self._m.append(m)
        self._m = tuple(self._m)
        SurfaceMapping.__init__(self, self._m[0].domain(), self._m[-1].codomain())

    def push_vector_forward(self, tangent_vector):
        r"""Applies the mapping to the provided vector."""
        for mapping in self._m:
            tangent_vector = mapping.push_vector_forward(tangent_vector)
        return tangent_vector

    def factors(self):
        r"""
        Return the factors of the map as a tuple ordered in the order of the application.
        So for example, a return of `(f, g, h)` would indicate the map $h \circ g \circ f$.
        """
        return self._m

    def is_invertible(self):
        r"""Return if this mapping is invertible."""
        for mapping in self._m:
            if not mapping.is_invertible():
                return False
        return True

    def is_locally_affine(self):
        r"""
        Return `True` if the map is known to be affine in local coordinates.

        We return True if all the maps in the composition are locally affine.
        """
        for mapping in self._m:
            if not mapping.is_locally_affine():
                return False
        return True

    def is_cellular(self):
        r"""
        Return `True` if the map is known to be cellular.

        That is, polygons defining the surface in the domain are mapped
        bijectively to polygons in the codomain.

        We return `True` if all the maps being composed are cellular.
        """
        for mapping in self._m:
            if not mapping.is_cellular():
                return False
        return True

    def is_orientation_preserving(self):
        r"""
        Return `True` if the mapping is orientation-preserving and
        `False` if not.

        By default we return `True`.
        """
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.

        Raises a NotImplementedError by default.
        """
        raise NotImplementedError('__invert__ has not been implemented.')


    def __invert__(self):
        r"""
        Return the inverse mapping.

        Raises a NotImplementedError if the map is not invertible.
        """
        if self.is_invertible():
            return SurfaceMappingComposition(*[~m for m in reversed(self._m)])
        else:
            raise NotImplementedError('This SurfaceMappingComposition is not invertible')

    def __str__(self):
        return f'Composition of {len(self._m)} surface mappings.'

    def __repr__(self):
        return f'Composition of {len(self._m)} surface mappings.'

class CellularMapping(SurfaceMapping):
    r"""Abstract class for a cellular mapping between surfaces.

    A cellular mapping is one that sends the polygons making up
    the domain of the surface to polygons making up the codomain
    of the surface.
    """
    def __init__(self, domain, codomain, ring = None, locally_affine = False):
        r"""
        Construct a CellularMapping. Note this class is abstract.

        Parameters are the domain surface, the codomain surface, the `ring` such
        that all maps between polygons are restrictions of affine maps defined
        over the given `ring`.
        """
        if ring is None:
            if domain.base_ring() == codomain.base_ring():
                self._ring = domain.base_ring()
            else:
                from sage.structure.element import get_coercion_model
                cm = get_coercion_model()
                self._ring = cm.common_parent(domain.base_ring(), codomain.base_ring())
        else:
            if not ring.has_coerce_map_from(domain.base_ring()):
                raise ValueError(
                    'The ring used must have a coersion from the base ring of the domain surface.'
                )
            self._ring = ring
        self._locally_affine = locally_affine
        self._codomain_tb = codomain.tangent_bundle(ring)
        SurfaceMapping.__init__(self, domain, codomain)

    # Implement the following two methods:
    def image_vertex(self, domain_label, domain_vertex):
        r"""
        The pair `(domain_label, domain_vertex)` determines a vertex
        of a polygon. This maps to a vertex of a polygon in the codomain,
        and the corresponding pair is returned.
        """
        raise NotImplementedError

    def affine_transformation(self, domain_label):
        r"""
        Return the affine transformation carrying the polygon
        with label `domain_label` to the image polygon.

        The returned element should be in the AffineGroup over
        the provided ring. See :meth:`~flatsurf.geometry.mapping.CellularMapping.affine_group`.
        """
        raise NotImplementedError

    def ring(self):
        return self._ring

    def affine_group(self):
        return AffineGroup(2, self._ring)

    def is_locally_affine(self):
        r"""
        Return `True` if the map is affine in local coordinates.

        Note that this requires that the map is affine at points in the
        interiors of edges.
        """
        return self._locally_affine

    def is_cellular(self):
        r"""
        Return `True` if the map is known to be cellular.

        That is, polygons defining the surface in the domain are mapped
        bijectively to polygons in the codomain.
        """
        return True

    def _test_vertices(self, **options):
        r"""
        Check that the map sends vertices to the vertices it claims it does.
        """
        tester = self._tester(**options)
        count = 0
        for l in self.domain().label_iterator():
            p = self.domain().polygon(l)
            ne = p.num_edges()
            g = self.affine_transformation(l)
            ll = self.image_vertex(l,0)[0]
            pp = self.codomain().polygon(ll)
            tester.assertEqual(ne, pp.num_edges(),
                f'A cellular map can not send a polygon with {ne} edges to one with {pp.num_edges()}! Domain label is {l}.')
            for v in range(ne):
                vv = self.image_vertex(l,v)[1]
                tester.assertEqual(g(p.vertex(v)), pp.vertex(vv),
                    f'The map is supposed to send vertex ({l},{v}) to ({ll}, {vv}), but does not.')
            count += 1
            if (not self.domain().is_finite()) and count >= 30:
                break

    def image_edge(self, domain_label, domain_edge):
        r"""
        Return the pair
        """
        if self.is_orientation_preserving():
            return self.image_vertex(domain_label, domain_edge)
        else:
            ne = self.domain().polygon(domain_label).num_edges()
            return self.image_vertex(domain_label, (domain_edge+1)%ne)

    def affine_transformation(self, domain_label):
        r"""
        Return the affine transformation carrying the polygon
        with label `domain_label` to the image polygon.

        The returned element should be in the AffineGroup over
        the provided ring. See :meth:`~flatsurf.geometry.mapping.CellularMapping.affine_group`.
        """
        raise NotImplementedError

    def push_vector_forward(self, tangent_vector, ring=None):
        r"""
        Applies the mapping to the provided vector.
        """
        domain_label = tangent_vector.polygon_label()
        codomain_label = self.image_vertex(domain_label,0)[0]
        at = self.affine_transformation(domain_label)
        if ring is None:
            return self._codomain_tb(
                codomain_label,
                at(tangent_vector.point()),
                at.matrix()[:2,:2]*tangent_vector.vector()
            )
        else:
            return self.codomain().tangent_bundle(ring)(
                codomain_label,
                at(tangent_vector.point()),
                at.matrix()[:2,:2]*tangent_vector.vector()
            )

    def _test_affine_along_edges(self, **options):
        r"""
        Check that the map is cellular: i.e., that it sends polygons
        to polygons.
        """
        tester = self._tester(**options)
        count = 0
        AG = self.affine_group()
        for l,e in self.domain().edge_iterator():
            g = self.affine_transformation(l)
            ll, ee = self.domain().opposite_edge(l,e)
            gg = self.affine_transformation(ll)
            t1 = self.domain().edge_transformation(l, e)
            image_l,image_e = self.image_edge(l, e)
            t2 = self.codomain().edge_transformation(image_l, image_e)
            tester.assertEqual(gg*AG(t1), AG(t2)*g,
                f'The mapping is not affine along edge ({l},{e}).')
            count += 1
            if (not self.domain().is_finite()) and count >= 30:
                break

class IdentityMapping(CellularMapping):
    r"""
    Construct an identity map between two "equal" surfaces.

    Surfaces should have the same labels, the same polygons, and the same
    gluings, but the surfaces may be defined over slightly different rings.
    The map will be defined over the common parent.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s1 = translation_surfaces.regular_octagon()
        sage: s2 = s1.copy(new_field=AA)
        sage: from flatsurf.geometry.mappings import IdentityMapping
        sage: m = IdentityMapping(s1, s2)
        sage: m.ring() == s2.base_ring()
        True
        sage: TestSuite(m).run()

        sage: m2 = IdentityMapping(s1, s1)
        sage: TestSuite(m2).run()
    """
    def __init__(self, domain, codomain = None, ring = None):
        if codomain is None:
            codomain = domain
        CellularMapping.__init__(self, domain, codomain, ring=ring, locally_affine = True)

    def _test_equality_of_domain_and_codomain(self, **options):
        r"""
        Check that the map is cellular: i.e., that it sends polygons
        to polygons.
        """
        tester = self._tester(**options)
        s1 = self.domain()
        s2 = self.codomain()
        if s1.base_ring() != s2.base_ring():
            r = self.ring()
            if s1 != r:
                s1 = s1.copy(new_field = r)
            if s2 != r:
                s2 = s2.copy(new_field = r)
        tester.assertEqual(s1, s2,
                f'The domain and codomain must be equal in the identity map up to a change of base_ring.')

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self._codomain.tangent_vector( \
            tangent_vector.polygon_label(), \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.
        """
        return self

    def image_vertex(self, domain_label, domain_vertex):
        r"""
        The pair `(domain_label, domain_vertex)` determines a vertex
        of a polygon. This maps to a vertex of a polygon in the codomain,
        and the corresponding pair is returned.
        """
        return (domain_label, domain_vertex)


    def affine_transformation(self, domain_label):
        r"""
        Return the affine transformation carrying the polygon
        with label `domain_label` to the image polygon.

        The returned element should be in the AffineGroup over
        the provided ring. See :meth:`~flatsurf.geometry.mapping.CellularMapping.affine_group`.
        """
        return self.affine_group().one()

class LCCSCanonicalizationMapping(CellularMapping):
    r'''
    This class represents the natural mapping from a surface to its
    :class:`~flatsurf.geometry.surface.LazyCellularCanonicalizeSurface`.

    Typically you would get this map from the
    :meth:`~flatsurf.geometry.surface.LazyCellularCanonicalizeSurface.canonicalization_mapping`
    and tests are in that function declaration.
    '''
    def __init__(self, lccs, codomain):
        from .surface import LazyCellularCanonicalizeSurface
        # lccs should be the The LazyCellularCanonicalizeSurface
        if not isinstance(lccs, LazyCellularCanonicalizeSurface):
            raise ValueError('Parameter `lccs` must be an instance of LazyCellularCanonicalizeSurface.')
        self._lccs = lccs # The LazyCellularCanonicalizeSurface
        CellularMapping.__init__(self, self._lccs._s, codomain, locally_affine = True)

    def image_vertex(self, domain_label, domain_vertex):
        r"""
        The pair `(domain_label, domain_vertex)` determines a vertex
        of a polygon. This maps to a vertex of a polygon in the codomain,
        and the corresponding pair is returned.
        """
        p = self._lccs._s.polygon(domain_label)
        il = self._lccs._found_labels[domain_label]
        zv = self._lccs._old_data[il][1] # Vertex in p mapping to the zero vertex
        ne = p.num_edges()
        return (il, (ne + domain_vertex - zv)%ne)

    def affine_transformation(self, domain_label):
        r"""
        Return the affine transformation carrying the polygon
        with label `domain_label` to the image polygon.

        The returned element should be in the AffineGroup over
        the provided ring. See :meth:`~flatsurf.geometry.mapping.CellularMapping.affine_group`.
        """
        i,v = self.image_vertex(domain_label, 0)
        zv = self.domain().polygon(domain_label).num_edges() - v # This vertex maps to 0.
        return self.affine_group()(self._lccs._matrices[i]) * \
            self.affine_group().translation(-self._domain.polygon(domain_label).vertex(zv))

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.
        """
        return self._lccs.inverse_mapping(self.codomain())

class LCCSCanonicalizationInverseMapping(CellularMapping):
    r'''
    This class represents the natural mapping from a surface to its
    :class:`~flatsurf.geometry.surface.LazyCellularCanonicalizeSurface`.

    Typically you would get this map from the
    :meth:`~flatsurf.geometry.surface.LazyCellularCanonicalizeSurface.canonicalization_mapping`
    and tests are in that function declaration.
    '''
    def __init__(self, lccs, domain):
        from .surface import LazyCellularCanonicalizeSurface
        # lccs should be the The LazyCellularCanonicalizeSurface
        if not isinstance(lccs, LazyCellularCanonicalizeSurface):
            raise ValueError('Parameter `lccs` must be an instance of LazyCellularCanonicalizeSurface.')
        self._lccs = lccs # The LazyCellularCanonicalizeSurface
        CellularMapping.__init__(self, domain, self._lccs._s, locally_affine = True)

    def image_vertex(self, domain_label, domain_vertex):
        r"""
        The pair `(domain_label, domain_vertex)` determines a vertex
        of a polygon. This maps to a vertex of a polygon in the codomain,
        and the corresponding pair is returned.
        """
        image_label, zv = self._lccs._old_data[domain_label]
        # zv is the label in the original surface (self._lccs._s) that
        # corresponds to the 0 vertex on the canonicalized surface.
        ne = self.domain().polygon(domain_label).num_edges()
        return (image_label, (domain_vertex+zv)%ne)

    def affine_transformation(self, domain_label):
        r"""
        Return the affine transformation carrying the polygon
        with label `domain_label` to the image polygon.

        The returned element should be in the AffineGroup over
        the provided ring. See :meth:`~flatsurf.geometry.mapping.CellularMapping.affine_group`.
        """
        image_label, zv = self._lccs._old_data[domain_label]
        return \
            self.affine_group().translation(self._codomain.polygon(image_label).vertex(zv)) * \
            self.affine_group()( ~ self._lccs._matrices[domain_label] )

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.
        """
        return self._lccs.mapping(self.domain())


class JoinPolygonsMapping(SurfaceMapping):
    r"""
    Return a SurfaceMapping joining two polygons together along the edge provided to the constructor.

    EXAMPLES::

        sage: from flatsurf.geometry.surface import Surface_list
        sage: from flatsurf.geometry.translation_surface import TranslationSurface
        sage: from flatsurf.geometry.polygon import ConvexPolygons
        sage: P = ConvexPolygons(QQ)
        sage: s0=Surface_list(base_ring=QQ)
        sage: s0.add_polygon(P([(1,0),(0,1),(-1,-1)])) # gets label=0
        0
        sage: s0.add_polygon(P([(-1,0),(0,-1),(1,1)])) # gets label=1
        1
        sage: s0.change_polygon_gluings(0,[(1,0),(1,1),(1,2)])
        sage: s0.set_immutable()
        sage: s=TranslationSurface(s0)
        sage: from flatsurf.geometry.mappings import *
        sage: m=JoinPolygonsMapping(s,0,2)
        sage: s2=m.codomain()
        sage: for label,polygon in s2.label_iterator(polygons=True):
        ....:     print("Polygon "+str(label)+" is "+str(polygon)+".")
        Polygon 0 is Polygon: (0, 0), (1, 0), (1, 1), (0, 1).
        sage: for label,edge in s2.edge_iterator():
        ....:     print(str((label,edge))+" is glued to "+str(s2.opposite_edge(label,edge))+".")
        (0, 0) is glued to (0, 2).
        (0, 1) is glued to (0, 3).
        (0, 2) is glued to (0, 0).
        (0, 3) is glued to (0, 1).
    """
    def __init__(self, s, p1, e1):
        r"""
        Join polygon with label p1 of s to polygon sharing edge e1.
        """
        if s.is_mutable():
            raise ValueError("Can only construct JoinPolygonsMapping for immutable surfaces.")

        ss2=s.copy(lazy=True,mutable=True)
        s2=ss2.underlying_surface()

        poly1=s.polygon(p1)
        p2,e2 = s.opposite_edge(p1,e1)
        poly2=s.polygon(p2)
        t=s.edge_transformation(p2, e2)
        dt=t.derivative()
        vs = []  # actually stores new edges...
        edge_map={} # Store the pairs for the old edges.
        for i in range(e1):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))
        ne=poly2.num_edges()
        for i in range(1,ne):
            ee=(e2+i)%ne
            edge_map[len(vs)]=(p2,ee)
            vs.append(dt * poly2.edge( ee ))
        for i in range(e1+1, poly1.num_edges()):
            edge_map[len(vs)]=(p1,i)
            vs.append(poly1.edge(i))

        inv_edge_map={}
        for key, value in iteritems(edge_map):
            inv_edge_map[value]=(p1,key)

        if s.base_label()==p2:
            # The polygon with the base label is being removed.
            s2.change_base_label(p1)

        s2.change_polygon(p1, ConvexPolygons(s.base_ring())(vs).translate(poly1.vertex(0)))

        for i in range(len(vs)):
            p3,e3 = edge_map[i]
            p4,e4 = s.opposite_edge(p3,e3)
            if p4 == p1 or p4 == p2:
                pp,ee = inv_edge_map[(p4,e4)]
                s2.change_edge_gluing(p1,i,pp,ee)
            else:
                s2.change_edge_gluing(p1,i,p4,e4)

        s2.set_immutable()

        self._saved_label=p1
        self._removed_label=p2
        self._remove_map = t
        self._remove_map_derivative = dt
        self._glued_edge=e1
        SurfaceMapping.__init__(self, s, ss2)

    def removed_label(self):
        r"""
        Return the label that was removed in the joining process.
        """
        return self._removed_label

    def glued_vertices(self):
        r"""
        Return the vertices of the newly glued polygon which bound the diagonal formed by the glue.
        """
        return (self._glued_edge,
                (self._glued_edge+self._domain.polygon(self._removed_label).num_edges()-1) %
                     self._codomain.polygon(self._saved_label).num_edges())

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._removed_label:
            return self._codomain.tangent_vector( \
                self._saved_label, \
                self._remove_map(tangent_vector.point()), \
                self._remove_map_derivative*tangent_vector.vector(), \
                ring = ring)
        else:
            return self._codomain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

#    def pull_vector_back(self,tangent_vector):
#        r"""
#        Applies the inverse of the mapping to the provided vector.
#        """
#        ring = tangent_vector.bundle().base_ring()
#        if tangent_vector.polygon_label() == self._saved_label:
#            p=tangent_vector.point()
#            v=self._domain.polygon(self._saved_label).vertex(self._glued_edge)
#            e=self._domain.polygon(self._saved_label).edge(self._glued_edge)
#            from flatsurf.geometry.polygon import wedge_product
#            wp = wedge_product(p-v,e)
#            if wp > 0:
#                # in polygon with the removed label
#                return self.domain().tangent_vector( \
#                    self._removed_label, \
#                    (~ self._remove_map)(tangent_vector.point()), \
#                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
#                    ring = ring)
#            if wp < 0:
#                # in polygon with the removed label
#                return self.domain().tangent_vector( \
#                    self._saved_label, \
#                    tangent_vector.point(), \
#                    tangent_vector.vector(), \
#                    ring = ring)
#            # Otherwise wp==0
#            w = tangent_vector.vector()
#            wp = wedge_product(w,e)
#            if wp > 0:
#                # in polygon with the removed label
#                return self.domain().tangent_vector( \
#                    self._removed_label, \
#                    (~ self._remove_map)(tangent_vector.point()), \
#                    (~ self._remove_map_derivative)*tangent_vector.vector(), \
#                    ring = ring)
#            return self.domain().tangent_vector( \
#                self._saved_label, \
#                tangent_vector.point(), \
#                tangent_vector.vector(), \
#                ring = ring)
#        else:
#            return self._domain.tangent_vector( \
#                tangent_vector.polygon_label(), \
#                tangent_vector.point(), \
#                tangent_vector.vector(), \
#                ring = ring)

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

    def __invert__(self):
        r"""
        Return the inverse mapping.

        Raises a NotImplementedError by default.
        """
        v1, v2 = self.glued_vertices()
        print((v1,v2))
        return SplitPolygonsMapping(self.codomain(), self._saved_label, v1, v2, new_label = self.removed_label())

    def is_locally_affine(self):
        r"""
        Return `True` if the map is affine in local coordinates.

        Note that this requires that the map is affine at points in the
        interiors of edges.
        """
        return True

    def is_cellular(self):
        r"""
        Return `True` if the map is cellular, i.e., polygons defining the
        surface in the domain are mapped bijectively to polygons in the
        codomain.
        """
        return True

    def is_orientation_preserving(self):
        r"""
        Return `True` if the mapping is orientation-preserving and
        `False` if not.

        By default we return `True`.
        """
        return True

class SplitPolygonsMapping(SurfaceMapping):
    r"""
    Class for cutting a polygon along a diagonal.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.mappings import SplitPolygonsMapping
        sage: m = SplitPolygonsMapping(s,0,0,2)
        sage: s2=m.codomain()
        sage: TestSuite(s2).run()
        sage: for pair in s2.label_iterator(polygons=True):
        ....:     print(pair)
        (0, Polygon: (0, 0), (1/2*a + 1, 1/2*a), (1/2*a + 1, 1/2*a + 1), (1, a + 1), (0, a + 1), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a))
        (ExtraLabel(...), Polygon: (0, 0), (-1/2*a - 1, -1/2*a), (-1/2*a, -1/2*a))
        sage: for glue in s2.edge_iterator(gluings=True):
        ....:     print(glue)
        ((0, 0), (ExtraLabel(...), 0))
        ((0, 1), (0, 5))
        ((0, 2), (0, 6))
        ((0, 3), (ExtraLabel(...), 1))
        ((0, 4), (ExtraLabel(...), 2))
        ((0, 5), (0, 1))
        ((0, 6), (0, 2))
        ((ExtraLabel(...), 0), (0, 0))
        ((ExtraLabel(...), 1), (0, 3))
        ((ExtraLabel(...), 2), (0, 4))
    """

    def __init__(self, s, p, v1, v2, new_label = None):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.

        Warning: We do not ensure that new_label is not already in the list of labels unless it is None (as by default).
        """
        if s.is_mutable():
            raise ValueError("The surface should be immutable.")

        poly=s.polygon(p)
        ne=poly.num_edges()
        if v1<0 or v2<0 or v1>=ne or v2>=ne:
            raise ValueError('Provided vertices out of bounds.')
        if abs(v1-v2)<=1 or abs(v1-v2)>=ne-1:
            raise ValueError('Provided diagonal is not a diagonal.')
        if v2<v1:
            temp=v1
            v1=v2
            v2=temp

        newvertices1=[poly.vertex(v2)-poly.vertex(v1)]
        for i in range(v2, v1+ne):
            newvertices1.append(poly.edge(i))
        newpoly1 = ConvexPolygons(s.base_ring())(newvertices1)

        newvertices2=[poly.vertex(v1)-poly.vertex(v2)]
        for i in range(v1,v2):
            newvertices2.append(poly.edge(i))
        newpoly2 = ConvexPolygons(s.base_ring())(newvertices2)

        ss2 = s.copy(mutable=True,lazy=True)
        s2 = ss2.underlying_surface()
        s2.change_polygon(p,newpoly1)
        new_label = s2.add_polygon(newpoly2, label=new_label)

        old_to_new_labels={}
        for i in range(ne):
            if i<v1:
                old_to_new_labels[i]=(p,i+ne-v2+1)
            elif i<v2:
                old_to_new_labels[i]=(new_label,i-v1+1)
            else: # i>=v2
                old_to_new_labels[i]=(p,i-v2+1)
        new_to_old_labels={}
        for i,pair in iteritems(old_to_new_labels):
            new_to_old_labels[pair]=i

        # This glues the split polygons together.
        s2.change_edge_gluing(p,0,new_label,0)
        for e in range(ne):
            ll,ee = old_to_new_labels[e]
            lll,eee = s.opposite_edge(p,e)
            if lll == p:
                gl,ge = old_to_new_labels[eee]
                s2.change_edge_gluing(ll,ee,gl,ge)
            else:
                s2.change_edge_gluing(ll,ee,lll,eee)

        s2.set_immutable()

        self._p=p
        self._v1=v1
        self._v2=v2
        self._new_label=new_label
        from flatsurf.geometry.similarity import SimilarityGroup
        TG = SimilarityGroup(s.base_ring())
        self._tp = TG(-s.polygon(p).vertex(v1))
        self._tnew_label = TG(-s.polygon(p).vertex(v2))
        SurfaceMapping.__init__(self, s, ss2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            point=tangent_vector.point()
            vertex1=self._domain.polygon(self._p).vertex(self._v1)
            vertex2=self._domain.polygon(self._p).vertex(self._v2)

            wp = wedge_product(vertex2-vertex1,point-vertex1)

            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector( \
                    self._p, \
                    self._tp(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)
            if wp < 0:
                # in new polygon 2
                return self.codomain().tangent_vector( \
                    self._new_label, \
                    self._tnew_label(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)

            # Otherwise wp==0
            w = tangent_vector.vector()
            wp = wedge_product(vertex2-vertex1,w)
            if wp > 0:
                # in new polygon 1
                return self.codomain().tangent_vector( \
                    self._p, \
                    self._tp(tangent_vector.point()), \
                    tangent_vector.vector(), \
                    ring = ring)
            # in new polygon 2
            return self.codomain().tangent_vector( \
                self._new_label, \
                self._tnew_label(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._codomain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)


    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        if tangent_vector.polygon_label() == self._p:
            return self._domain.tangent_vector( \
                self._p, \
                (~ self._tp)(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        elif tangent_vector.polygon_label() == self._new_label:
            return self._domain.tangent_vector( \
                self._p, \
                (~ self._tnew_label)(tangent_vector.point()), \
                tangent_vector.vector(), \
                ring = ring)
        else:
            # Not in a polygon that was changed. Just copy the data.
            return self._domain.tangent_vector( \
                tangent_vector.polygon_label(), \
                tangent_vector.point(), \
                tangent_vector.vector(), \
                ring = ring)

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

def subdivide_a_polygon(s):
    r"""
    Return a SurfaceMapping which cuts one polygon along a diagonal or None if the surface is triangulated.
    """
    from flatsurf.geometry.polygon import wedge_product
    for l,poly in s.label_iterator(polygons=True):
        n = poly.num_edges()
        if n>3:
            for i in range(n):
                e1=poly.edge(i)
                e2=poly.edge((i+1)%n)
                if wedge_product(e1,e2) != 0:
                    return SplitPolygonsMapping(s,l,i, (i+2)%n)
            raise ValueError("Unable to triangulate polygon with label "+str(l)+\
                ": "+str(poly))
    return None


def triangulation_mapping(s):
    r"""Return a  SurfaceMapping triangulating the provided surface.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s=translation_surfaces.veech_2n_gon(4)
        sage: from flatsurf.geometry.mappings import *
        sage: m=triangulation_mapping(s)
        sage: s2=m.codomain()
        sage: TestSuite(s2).run()
        sage: for label,polygon in s2.label_iterator(polygons=True):
        ....:     print(str(polygon))
        Polygon: (0, 0), (-1/2*a, 1/2*a + 1), (-1/2*a, 1/2*a)
        Polygon: (0, 0), (1/2*a, -1/2*a - 1), (1/2*a, 1/2*a)
        Polygon: (0, 0), (-1/2*a - 1, -1/2*a - 1), (0, -1)
        Polygon: (0, 0), (-1, -a - 1), (1/2*a, -1/2*a)
        Polygon: (0, 0), (0, -a - 1), (1, 0)
        Polygon: (0, 0), (-1/2*a - 1, -1/2*a), (-1/2*a, -1/2*a)
    """
    assert(s.is_finite())
    m=subdivide_a_polygon(s)
    if m is None:
        return None
    s1=m.codomain()
    while True:
        m2=subdivide_a_polygon(s1)
        if m2 is None:
            return m
        s1=m2.codomain()
        m=SurfaceMappingComposition(m,m2)
    return m

def flip_edge_mapping(s,p1,e1):
    r"""
    Return a mapping whose domain is s which flips the provided edge.
    """
    m1=JoinPolygonsMapping(s,p1,e1)
    v1,v2=m1.glued_vertices()
    removed_label = m1.removed_label()
    m2=SplitPolygonsMapping(m1.codomain(), p1, (v1+1)%4, (v1+3)%4, new_label = removed_label)
    return SurfaceMappingComposition(m1,m2)

def one_delaunay_flip_mapping(s):
    r"""
    Returns one delaunay flip, or none if no flips are needed.
    """
    for p,poly in s.label_iterator(polygons=True):
        for e in range(poly.num_edges()):
            if s._edge_needs_flip(p,e):
                return flip_edge_mapping(s,p,e)
    return None

def delaunay_triangulation_mapping(s):
    r"""
    Returns a mapping to a Delaunay triangulation or None if the surface already is Delaunay triangulated.
    """
    assert(s.is_finite())
    m=triangulation_mapping(s)
    if m is None:
        s1=s
    else:
        s1=m.codomain()
    m1=one_delaunay_flip_mapping(s1)
    if m1 is None:
        return m
    if m is None:
        m=m1
    else:
        m=SurfaceMappingComposition(m,m1)
    s1=m1.codomain()
    while True:
        m1=one_delaunay_flip_mapping(s1)
        if m1 is None:
            return m
        s1=m1.codomain()
        m=SurfaceMappingComposition(m,m1)

def delaunay_decomposition_mapping(s):
    r"""
    Returns a mapping to a Delaunay decomposition or possibly None if the surface already is Delaunay.
    """
    m=delaunay_triangulation_mapping(s)
    if m is None:
        s1=s
    else:
        s1=m.codomain()
    edge_vectors=[]
    lc = s._label_comparator()
    for p,poly in s1.label_iterator(polygons=True):
        for e in range(poly.num_edges()):
            pp,ee=s1.opposite_edge(p,e)
            if (lc.lt(p,pp) or (p==pp and e<ee)) and s1._edge_needs_join(p,e):
                edge_vectors.append( s1.tangent_vector(p,poly.vertex(e),poly.edge(e)) )
    if len(edge_vectors)>0:
        ev=edge_vectors.pop()
        p,e=ev.edge_pointing_along()
        m1=JoinPolygonsMapping(s1,p,e)
        s2=m1.codomain()
        while len(edge_vectors)>0:
            ev=edge_vectors.pop()
            ev2=m1.push_vector_forward(ev)
            p,e=ev2.edge_pointing_along()
            mtemp=JoinPolygonsMapping(s2,p,e)
            m1=SurfaceMappingComposition(m1,mtemp)
            s2=m1.codomain()
        if m is None:
            return m1
        else:
            return SurfaceMappingComposition(m,m1)
    return m

def canonical_first_vertex(polygon):
    r"""
    Return the index of the vertex with smallest y-coordinate.
    If two vertices have the same y-coordinate, then the one with least x-coordinate is returned.
    """
    best=0
    best_pt=polygon.vertex(best)
    for v in range(1,polygon.num_edges()):
        pt=polygon.vertex(v)
        if pt[1]<best_pt[1]:
            best=v
            best_pt=pt
    if best==0:
        if pt[1]==best_pt[1]:
            return v
    return best

class CanonicalizePolygonsMapping(SurfaceMapping):
    r"""
    This is a mapping to a surface with the polygon vertices canonically determined.
    A canonical labeling is when the canonocal_first_vertex is the zero vertex.
    """
    def __init__(self, s):
        r"""
        Split the polygon with label p of surface s along the diagonal joining vertex v1 to vertex v2.
        """
        if not s.is_finite():
            raise ValueError("Currently only works with finite surfaces.")
        ring=s.base_ring()
        from flatsurf.geometry.similarity import SimilarityGroup
        T = SimilarityGroup(ring)
        P = ConvexPolygons(ring)
        cv = {} # dictionary for canonical vertices
        translations={} # translations bringing the canonical vertex to the origin.
        s2 = Surface_dict(base_ring=ring)
        for l,polygon in s.label_iterator(polygons=True):
            cv[l]=cvcur=canonical_first_vertex(polygon)
            newedges=[]
            for i in range(polygon.num_edges()):
                newedges.append(polygon.edge( (i+cvcur) % polygon.num_edges() ))
            s2.add_polygon(P(newedges), label=l)
            translations[l]=T( -polygon.vertex(cvcur) )
        for l1,polygon in s.label_iterator(polygons=True):
            for e1 in range(polygon.num_edges()):
                l2,e2=s.opposite_edge(l1,e1)
                ee1= (e1-cv[l1]+polygon.num_edges())%polygon.num_edges()
                polygon2=s.polygon(l2)
                ee2= (e2-cv[l2]+polygon2.num_edges())%polygon2.num_edges()
                # newgluing.append( ( (l1,ee1),(l2,ee2) ) )
                s2.change_edge_gluing(l1,ee1,l2,ee2)
        s2.change_base_label(s.base_label())
        s2.set_immutable()
        ss2=s.__class__(s2)

        self._cv=cv
        self._translations=translations

        SurfaceMapping.__init__(self, s, ss2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        l=tangent_vector.polygon_label()
        return self.codomain().tangent_vector(l, \
            self._translations[l](tangent_vector.point()), \
            tangent_vector.vector(), \
            ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        l=tangent_vector.polygon_label()
        return self.domain().tangent_vector(l, \
            (~self._translations[l])(tangent_vector.point()), \
            tangent_vector.vector(), \
            ring = ring)

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True

class ReindexMapping(SurfaceMapping):
    r"""
    Apply a dictionary to relabel the polygons.
    """
    def __init__(self,s,relabler,new_base_label=None):
        r"""
        The parameters should be a surface and a dictionary which takes as input a label and produces a new label.
        """
        if not s.is_finite():
            raise ValueError("Currently only works with finite surfaces.""")
        f = {} # map for labels going forward.
        b = {} # map for labels going backward.
        for l in s.label_iterator():
            if l in relabler:
                l2=relabler[l]
                f[l]=l2
                if l2 in b:
                    raise ValueError("Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change.")
                b[l2]=l
            else:
                # If no key then don't change the label
                f[l]=l
                if l in b:
                    raise ValueError("Provided dictionary has two keys mapping to the same value. Or you are mapping to a label you didn't change.")
                b[l]=l

        self._f=f
        self._b=b

        if new_base_label==None:
            if s.base_label() in f:
                new_base_label = f[s.base_label()]
            else:
                new_base_label = s.base_label()
        s2=s.copy(mutable=True,lazy=True)
        s2.relabel(relabler, in_place=True)
        s2.underlying_surface().change_base_label(new_base_label)

        SurfaceMapping.__init__(self, s, s2)

    def push_vector_forward(self,tangent_vector):
        r"""Applies the mapping to the provided vector."""
        # There is no change- we just move it to the new surface.
        ring = tangent_vector.bundle().base_ring()
        return self.codomain().tangent_vector( \
            self._f[tangent_vector.polygon_label()], \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

    def pull_vector_back(self,tangent_vector):
        r"""Applies the pullback mapping to the provided vector."""
        ring = tangent_vector.bundle().base_ring()
        return self.domain().tangent_vector( \
            self._b[tangent_vector.polygon_label()], \
            tangent_vector.point(), \
            tangent_vector.vector(), \
            ring = ring)

    def is_invertible(self):
        r"""Return true to indicate this mapping is invertible."""
        return True


def my_sgn(val):
    if val>0:
        return 1
    elif val<0:
        return -1
    else:
        return 0

def polygon_compare(poly1,poly2):
    r"""
    Compare two polygons first by area, then by number of sides,
    then by lexigraphical ording on edge vectors."""
    # This should not be used is broken!!
    #from sage.functions.generalized import sgn
    res = my_sgn(-poly1.area()+poly2.area())
    if res!=0:
        return res
    res = my_sgn(poly1.num_edges()-poly2.num_edges())
    if res!=0:
        return res
    ne=poly1.num_edges()
    for i in range(0,ne-1):
        edge_diff = poly1.edge(i) - poly2.edge(i)
        res = my_sgn(edge_diff[0])
        if res!=0:
            return res
        res = my_sgn(edge_diff[1])
        if res!=0:
            return res
    return 0

def translation_surface_cmp(s1, s2):
    r"""
    Compare two finite surfaces.
    The surfaces will be considered equal if and only if there is a translation automorphism
    respecting the polygons and the base_labels.
    """
    if not s1.is_finite() or not s2.is_finite():
        raise NotImplementedError
    lw1=s1.walker()
    lw2=s2.walker()
    try:
        from itertools import zip_longest
    except ImportError:
        from itertools import izip_longest as zip_longest
    for p1,p2 in zip_longest(lw1.polygon_iterator(), lw2.polygon_iterator()):
        if p1 is None:
            # s2 has more polygons
            return -1
        if p2 is None:
            # s1 has more polygons
            return 1
        ret = polygon_compare(p1,p2)
        if ret != 0:
            return ret
    # Polygons are identical. Compare edge gluings.
    for pair1,pair2 in zip_longest(lw1.edge_iterator(), lw2.edge_iterator()):
        l1,e1 = s1.opposite_edge(pair1)
        l2,e2 = s2.opposite_edge(pair2)
        num1 = lw1.label_to_number(l1)
        num2 = lw2.label_to_number(l2)
        ret = (num1 > num2) - (num1 < num2)
        if ret!=0:
            return ret
        ret = (e1 > e2) - (e1 < e2)
        if ret!=0:
            return ret
    return 0

def canonicalize_translation_surface_mapping(s):
    r"""
    Return the translation surface in a canonical form.

    EXAMPLES::

        sage: from flatsurf import *
        sage: s=translation_surfaces.octagon_and_squares().canonicalize()
        sage: TestSuite(s).run()
        sage: a = s.base_ring().gen()  # a is the square root of 2.

        sage: from flatsurf.geometry.mappings import *
        sage: mat=Matrix([[1,2+a],[0,1]])
        sage: from flatsurf.geometry.half_dilation_surface import GL2RMapping
        sage: m1=GL2RMapping(s,mat)
        sage: m2=canonicalize_translation_surface_mapping(m1.codomain())
        sage: m=m2*m1
        sage: translation_surface_cmp(m.domain(),m.codomain())==0
        True
        sage: TestSuite(m.codomain()).run()
        sage: s=m.domain()
        sage: v=s.tangent_vector(0,(0,0),(1,1))
        sage: w=m.push_vector_forward(v)
        sage: print(w)
        SimilaritySurfaceTangentVector in polygon 0 based at (0, 0) with vector (a + 3, 1)
    """
    from flatsurf.geometry.translation_surface import TranslationSurface
    if not s.is_finite():
        raise NotImplementedError
    if not isinstance(s,TranslationSurface):
        raise ValueError("Only defined for TranslationSurfaces")
    m1=delaunay_decomposition_mapping(s)
    if m1 is None:
        s2=s
    else:
        s2=m1.codomain()
    m2=CanonicalizePolygonsMapping(s2)
    if m1 is None:
        m=m2
    else:
        m=SurfaceMappingComposition(m1,m2)
    s2=m.codomain()

    s2copy=s2.copy(mutable=True)
    ss=s2.copy(mutable=True)
    labels={label for label in s2.label_iterator()}
    labels.remove(s2.base_label())
    for label in labels:
        ss.underlying_surface().change_base_label(label)
        if ss.cmp(s2copy)>0:
            s2copy.underlying_surface().change_base_label(label)
    # We now have the base_label correct.
    # We will use the label walker to generate the canonical labeling of polygons.
    w=s2copy.walker()
    w.find_all_labels()

    m3=ReindexMapping(s2,w.label_dictionary(),0)
    return SurfaceMappingComposition(m,m3)
