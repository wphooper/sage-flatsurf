from __future__ import absolute_import, print_function, division
from six.moves import range, map, filter, zip

from sage.misc.cachefunc import cached_method

from sage.structure.element import MultiplicativeGroupElement, parent
from sage.structure.unique_representation import UniqueRepresentation

from sage.categories.groups import Groups

from sage.all import Fields

from sage.modules.free_module_element import vector
from sage.groups.group import Group
from sage.rings.integer import Integer
from sage.rings.integer_ring import ZZ
from sage.modules.free_module_element import FreeModuleElement

from sage.env import SAGE_VERSION
if SAGE_VERSION >= '8.2':
    from sage.structure.element import is_Matrix
else:
    from sage.matrix.matrix import is_Matrix

from flatsurf.geometry.polygon import ConvexPolygon, ConvexPolygons

ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_m1 = -ZZ_1

class Similarity(MultiplicativeGroupElement):
    r"""
    Class for a similarity of the plane with possible reflection.

    Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t) if sign=1,
    and (ax+by+s,bx-ay+t) if sign=-1
    """
    def __init__(self, p, a, b, s, t, sign):
        r"""
        Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t) if sign=1,
        and (ax+by+s,bx-ay+t) if sign=-1
        """
        if p is None:
            raise ValueError("The parent must be provided")
        field = p._field
        if parent(a) is not field or \
           parent(b) is not field or \
           parent(s) is not field or \
           parent(t) is not field:
               raise ValueError("wrong parent for a,b,s or t")
        self._a = a
        self._b = b
        self._s = s
        self._t = t
        if parent(sign) is not ZZ or not sign.is_unit():
            raise ValueError("sign must be either 1 or -1.")
        self._sign = sign
        MultiplicativeGroupElement.__init__(self, p)

    def sign(self):
        return self._sign

    def is_translation(self):
        r"""
        Return whether this element is a translation.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_translation()
            True
            sage: S((1,0,3,-1/2)).is_translation()
            True
            sage: S((0,1,0,0)).is_translation()
            False
        """
        return self._sign.is_one() and self._a.is_one() and self._b.is_zero()

    def is_half_translation(self):
        r"""
        Return whether this element is a half translation.

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_half_translation()
            True
            sage: S((-1, 0, 0, 2)).is_half_translation()
            True
            sage: S((0,1,0,0)).is_half_translation()
            False
        """
        return self._sign.is_one() and (self._a.is_one() or ((-self._a).is_one())) and self._b.is_zero()

    def is_orientable(self):
        return self._sign.is_one()

    def is_rotation(self):
        r"""
        Check whether this element is a rotation

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)).is_rotation()
            False
            sage: S((0,-1,0,0)).is_rotation()
            True
            sage: S.one().is_rotation()
            True
        """
        return self.is_one() or (self.det().is_one() and not self.is_translation())

    def is_isometry(self):
        r"""
        Check whether this element is an isometry

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S.one().is_isometry()
            True
            sage: S((0,1,0,0)).is_isometry()
            True
            sage: S((0,1,0,0,-1)).is_isometry()
            True
            sage: S((1,1,0,0)).is_isometry()
            False
            sage: S((3,-1/2)).is_isometry()
            True
        """
        det = self.det()
        return det.is_one() or (-det).is_one()

    def det(self):
        r"""
        Return the determinant of this element
        """
        return self._sign * (self._a*self._a + self._b*self._b)

    def _mul_(left, right):
        r"""
        Composition

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,2)) * S((3,-5)) == S((4,-3))
            True

            sage: from itertools import product
            sage: a1 = S((0,2,0,0,1))
            sage: a2 = S((1,0,0,0,-1))
            sage: a3 = S((1,1,0,0))
            sage: a4 = S((1,0,-1,1))
            sage: a5 = S((2,-1,3/5,2/3,-1))
            sage: for g1,g2,g3 in product([a1,a2,a3,a4,a5], repeat=3):
            ....:     assert g1.matrix()*g2.matrix() == (g1*g2).matrix()
            ....:     assert (g1*g2).matrix()*g3.matrix() == (g1*g2*g3).matrix()
        """
        a = left._a * right._a - left._sign * left._b * right._b
        b = left._b * right._a + left._sign * left._a * right._b
        s = left._a * right._s - left._sign * left._b * right._t + left._s
        t = left._b * right._s + left._sign * left._a * right._t + left._t
        sign = left._sign * right._sign
        P = left.parent()
        return P.element_class(P, a, b, s, t, sign)

    def __invert__(self):
        r"""
        Invert a similarity.

        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: from itertools import product
            sage: for a in [S((0,2,0,0,1)), S((1,0,0,0,-1)), S((1,1,0,0)),
            ....:           S((1,0,-1,1)), S((2,-1,3/5,2/3,-1))]:
            ....:     assert (a*~a).is_one() and (~a*a).is_one()
        """
        P = self.parent()
        sign = self._sign
        det = self.det()
        a = sign*self._a/det
        b = -self._b/det
        return P.element_class(P,a,b,
            -a*self._s + sign*b*self._t,
            -b*self._s - sign*a*self._t,
            sign)

    def _div_(self, s):
        return self._mul_(s.__invert__())

    def __hash__(self):
        return 73*hash(self._a)-19*hash(self._b)+13*hash(self._s)+53*hash(self._t)+67*hash(self._sign)

    def __call__(self,w, field = None):
        r"""
        Return the image of ``w`` under the similarity. Here ``w`` may be a ConvexPolygon or a vector
        (or something that can be indexed in the same way as a vector). If a field is provided,
        the objects returned will be defined over this field.

        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(AA)
            sage: a = S((1,-1,AA(2).sqrt(),0))
            sage: a((1,2))
            (4.414213562373095?, 1)
            sage: a.matrix()*vector((1,2,1))
            (4.414213562373095?, 1, 1)

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SG = SimilarityGroup(QQ)
            sage: from flatsurf import ConvexPolygons
            sage: P = ConvexPolygons(QQ)
            sage: p = P.an_element()
            sage: p
            Polygon: (0, 0), (1, 0), (1, 1), (0, 1)
            sage: g = SG.an_element()**2
            sage: g
            (x, y) |-> (25*x + 4, 25*y + 10)
            sage: g(p)
            Polygon: (4, 10), (29, 10), (29, 35), (4, 35)
            sage: g(p, field=AA).parent()
            ConvexPolygons(Algebraic Real Field)
        """
        if field is not None:
            if not field in Fields():
                raise TypeError("field must be a field")
        if isinstance(w,ConvexPolygon):
            if field is None:
                P = ConvexPolygons(self.parent().base_field())
            else:
                P = ConvexPolygons(field)
            try:
                return P(vertices=[self(v) for v in w.vertices()])
            except ValueError as e:
                if not self._sign.is_one():
                    raise ValueError("Similarity must be orientation preserving.")
                else:
                    # Not sure why this would happen:
                    raise e
        if field is None:
            if self._sign.is_one():
                return vector([self._a*w[0]-self._b*w[1]+self._s, self._b*w[0]+self._a*w[1]+self._t])
            else:
                return vector([self._a*w[0]+self._b*w[1]+self._s, self._b*w[0]-self._a*w[1]+self._t])
        else:
            if self._sign.is_one():
                return vector(field, [self._a*w[0]-self._b*w[1]+self._s, self._b*w[0]+self._a*w[1]+self._t])
            else:
                return vector(field, [self._a*w[0]+self._b*w[1]+self._s, self._b*w[0]-self._a*w[1]+self._t])

    def _repr_(self):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S.one()
            (x, y) |-> (x, y)
            sage: S((1,-2/3))
            (x, y) |-> (x + 1, y - 2/3)
            sage: S((-1,0,2/3,3))
            (x, y) |-> (-x + 2/3, -y + 3)
            sage: S((-1,0,2/3,3,-1))
            (x, y) |-> (-x + 2/3, y + 3)
        """
        R = self.parent()._field['x','y']
        x,y = R.gens()
        return "(x, y) |-> ({}, {})".format(
                    self._a*x - self._sign*self._b*y + self._s,
                    self._b*x + self._sign*self._a*y + self._t)

    def __eq__(self, other):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,0)) == S((1,0))
            True
            sage: S((1,0)) == S((0,1))
            False
            sage: S((1,0,0,0)) == S((0,1,0,0))
            False
            sage: S((1,0,0,0,1)) == S((1,0,0,0,-1))
            False
        """
        if other is None:
            return False
        if type(other)==int:
            return False
        if self.parent() != other.parent():
            return False
        return self._a == other._a and \
               self._b == other._b and \
               self._s == other._s and \
               self._t == other._t and \
               self._sign == other._sign

    def __ne__(self, other):
        return not (self == other)

    def matrix(self):
        r"""
        Return the 3x3 matrix representative of this element

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)

            sage: S((1,-2/3,1,1,-1)).matrix()
            [   1 -2/3    1]
            [-2/3   -1    1]
            [   0    0    1]
        """
        P = self.parent()
        M = P._matrix_space_3x3()
        z = P._field.zero()
        o = P._field.one()
        return M(
            [self._a, -self._sign*self._b, self._s,
             self._b, +self._sign*self._a, self._t,
            z, z, o])

    def derivative(self):
        r"""
        Return the 2x2 matrix corresponding to the derivative of this element

        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)

            sage: S((1,-2/3,1,1,-1)).derivative()
            [   1 -2/3]
            [-2/3   -1]
        """
        M = self.parent()._matrix_space_2x2()
        return M([self._a, -self._sign*self._b, self._b,  self._sign*self._a])

    # OLD AND DEPRECATED

    def a(self):
        from sage.misc.superseded import deprecation
        deprecation(42, "Do not use .a()")
        return self._a

    def b(self):
        from sage.misc.superseded import deprecation
        deprecation(42, "Do not use .b()")
        return self._b

    def s(self):
        from sage.misc.superseded import deprecation
        deprecation(42, "Do not use .s()")
        return self._s

    def t(self):
        from sage.misc.superseded import deprecation
        deprecation(42, "Do not use .t()")
        return self._t

class SimilarityGroup(UniqueRepresentation, Group):
    r"""
    The group of possibly orientation reversing similarities in the plane.

    This is the group generated by rotations, translations and dilations.
    """
    Element = Similarity

    def __init__(self, base_field):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: TestSuite(SimilarityGroup(QQ)).run()
            sage: TestSuite(SimilarityGroup(AA)).run()
        """
        if not base_field in Fields():
            raise TypeError("base_field must be a field")
        self._field = base_field
        # The vector space of vectors
        Group.__init__(self, category=Groups().Infinite())

    @cached_method
    def _matrix_space_2x2(self):
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(self._field, 2)

    @cached_method
    def _matrix_space_3x3(self):
        from sage.matrix.matrix_space import MatrixSpace
        return MatrixSpace(self._field, 3)

    @cached_method
    def _vector_space(self):
        from sage.modules.free_module import VectorSpace
        return VectorSpace(self._field, 2)

    def _element_constructor_(self, *args, **kwds):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: S = SimilarityGroup(QQ)
            sage: S((1,1))  # translation
            (x, y) |-> (x + 1, y + 1)

            sage: V = QQ^2
            sage: S(V((1,-1)))
            (x, y) |-> (x + 1, y - 1)

            sage: S(vector((1,1)))
            (x, y) |-> (x + 1, y + 1)
        """
        if len(args) == 1:
            x = args[0]
        else:
            x = args

        a = self._field.one()
        b = s = t = self._field.zero()
        sign = ZZ_1

        # TODO: 2x2 and 3x3 matrix input

        if isinstance(x, (tuple,list)):
            if len(x) == 2:
                s,t = map(self._field, x)
            elif len(x) == 4:
                a,b,s,t = map(self._field, x)
            elif len(x) == 5:
                a,b,s,t = map(self._field, x[:4])
                sign = ZZ(x[4])
            else:
                raise ValueError("can not construct a similarity from a list of length {}".format(len(x)))
        elif is_Matrix(x):
            #   a -sb
            #   b sa
            if x.nrows() == x.ncols() == 2:
                a,c,b,d = x.list()
                if a == d and b == -c:
                    sign = ZZ_1
                elif a == -d and b == c:
                    sign = ZZ_m1
                else:
                    raise ValueError("not a similarity matrix")
            elif x.nrows() == x.ncols() == 3:
                raise NotImplementedError
            else:
                raise ValueError("invalid dimension for matrix input")
        elif isinstance(x, FreeModuleElement):
            if len(x) == 2:
                if x.base_ring() is self._field:
                    s,t = x
                else:
                    s,t = map(self._field, x)
            else:
                raise ValueError("invalid dimension for vector input")
        else:
            p = parent(x)
            if self._field.has_coerce_map_from(p):
                a = self._field(x)
            else:
                raise ValueError

        if (a*a + b*b).is_zero():
            raise ValueError("not invertible")

        return self.element_class(self, a, b, s, t, sign)

    def _coerce_map_from_(self, S):
        if self._field.has_coerce_map_from(S):
            return True
        if isinstance(S, SimilarityGroup):
            return self._field.has_coerce_map_from(S._field)

    def _repr_(self):
        r"""
        TESTS::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SimilarityGroup(QQ)
            Similarity group over Rational Field
        """
        return "Similarity group over {}".format(self._field)

    def one(self):
        r"""
        EXAMPLES::

            sage: from flatsurf.geometry.similarity import SimilarityGroup
            sage: SimilarityGroup(QQ).one()
            (x, y) |-> (x, y)
            sage: SimilarityGroup(QQ).one().is_one()
            True
        """
        return self.element_class(self,
                self._field.one(),  # a
                self._field.zero(), # b
                self._field.zero(), # s
                self._field.zero(), # t
                ZZ_1)               # sign

    def an_element(self):
        return self(3, 4, 2, -1, -1)

    def is_abelian(self):
        return False

    def base_field(self):
        return self._field
