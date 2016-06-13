from sage.matrix.constructor import matrix
from sage.matrix.matrix_space import MatrixSpace
from sage.groups.group import Group
from sage.categories.groups import Groups
from sage.structure.element import MultiplicativeGroupElement, parent
from sage.modules.free_module import VectorSpace
from sage.modules.free_module_element import vector
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.integer import Integer

from sage.categories.fields import Fields
_Fields = Fields()

from flatsurf.geometry.translation import TranslationGroup

ZZ_0 = Integer(0)
ZZ_1 = Integer(1)
ZZ_2 = Integer(2)
ZZ_3 = Integer(3)
ZZ_4 = Integer(4)


class Similarity(MultiplicativeGroupElement):
    r"""Class for a similarity of the plane."""
    
    def __init__(self, parent, a, b, s, t):
        r'''Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t).'''
        if parent is None:
            raise ValueError("The parent must be provided")
        self._a=a
        self._b=b
        self._s=s
        self._t=t
        self._parent=parent
        MultiplicativeGroupElement.__init__(self,parent)

    def _mul_(self,s):
        r'''Compose two similarities.'''
        C = self.__class__
        return C(self._parent,
            self._a*s._a-self._b*s._b, 
            self._b*s._a+self._a*s._b,
            self._a*s._s-self._b*s._t+self._s,
            self._b*s._s+self._a*s._t+self._t)

    def __invert__(self):
        r'''Invert a similarity.'''
        det=self._a*self._a+self._b*self._b
        a=self._a/det
        b=-self._b/det
        C = self.__class__
        return C(self._parent,
            a,
            b,
            -a*self._s+b*self._t,
            -b*self._s-a*self._t)

    def _div_(self,s):
        return self._mul_(s.__invert__())

    def __hash__(self):
        return 73*hash(self._a)-19*hash(self._b)+13*hash(self._s)+53*hash(self._t)

    def __call__(self,w):
        r'''Return m*w+v.'''
        return vector([self._a*w[0]-self._b*w[1]+self._s, self._b*w[0]+self._a*w[1]+self._t])

    def _repr_(self):
        return "Similarity (x,y) mapsto ("+str(self._a)+"*x-"+\
            str(self._b)+"*y+"+str(self._s)+", "+\
            str(self._b)+"*x+"+str(self._a)+"*y+"+str(self._t)+")"

    def _cmp_(self, other):
        x=cmp(self._a,other._a)
        if x!=0:
            return x
        x=cmp(self._b,other._b)
        if x!=0:
            return x
        x=cmp(self._s,other._s)
        if x!=0:
            return x
        return cmp(self._t,other._t)

    __cmp__=_cmp_

    # For pickling:
    #def __reduce__(self):
    #    return self.__class__, (self._parent, self._a, self._b, self._s, self._t)

    def matrix(self):
        return matrix(self._parent._f,[
            [self._a, -self._b, self._s],
            [self._b,  self._a, self._t],
            [self._parent._f.zero(), self._parent._f.zero(), self._parent._f.one()]])

    def a(self):
        return self._a

    def b(self):
        return self._b

    def s(self):
        return self._s

    def t(self):
        return self._t
        
    def sign(self):
        return 1

    def derivative(self):
        r"""Return the 2x2 matrix corresponding to the derivative of the similarity of the plane."""
        return matrix(self._parent._f,[
            [self._a, -self._b],
            [self._b,  self._a]])

    def det(self):
        r"""
        Return the determinant of the derivative of the map.
        """
        return self._a**2+self._b**2

class SimilarityGroup(UniqueRepresentation,Group):
    r'''Group representing all similarities in the plane.
    This is the group generated by rotations, translations and dilations.
    '''

    Element = Similarity

    def _element_constructor_(self, *args, **kwds):
        if len(args)!=1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        p=parent(x)
        if self._f.has_coerce_map_from(p):
            return self.element_class( self,self._f(x), self._f.zero(), self._f.zero(), self._f.zero())
        if isinstance(p, SimilarityGroup):
            return self.element_class(self, x.a(), x.b(), x.s(), x.t())
        if isinstance(p, TranslationGroup):
            return self.element_class( self,self._f.one(), self._f.zero(), x.s(), x.t() )
        return self.element_class(self, x, **kwds)

    def _coerce_map_from_(self, S):
        if self._f.has_coerce_map_from(S):
            return True
        if isinstance(S, SimilarityGroup):
            return self._f.has_coerce_map_from(S._f)
        if isinstance(S, TranslationGroup):
            return self._f.has_coerce_map_from(S.base_field())
           
    def __init__(self, base_field):
        self._f=base_field
        # The vector space of vectors 
        self._vs = VectorSpace(self._f,2)
        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        return "SimilarityGroup over field "+str(self._f)

    def one(self):
        return self.element_class(self,self._f.one(),self._f.zero(),self._f.zero(),self._f.zero())

    def an_element(self):
        return self.element_class(self,self._f(ZZ_3),self._f(ZZ_4),self._f(ZZ_2),self._f(-ZZ_1))

    def is_abelian(self):
        return False

    def gens(self):
        pairs=[
            (self._f.one(),self._f.zero()),
            (self._f(ZZ_2),self._f.zero()),
            (self._f.zero(),self._f(ZZ_2)),
            (self._f(ZZ_3),self._f(ZZ_4))]
        l=[]
        for p in pairs:
            for v in self._vs.gens():
                l.append(self.element_class(self,p[0],p[1],v[0],v[1]))
        return l
    
    # For pickling:
    #def __reduce__(self):
    #    return self.__class__, (self._f,)
        
    #def _cmp_(self, other):
    #    return self._f == other._f

    #__cmp__=_cmp_

    def base_field(self):
        return self._f

class SimilarityReflection(MultiplicativeGroupElement):
    r"""Class for a similarity of the plane with possible reflection."""
    
    def __init__(self, parent, a, b, s, t, sign):
        r'''Construct the similarity (x,y) mapsto (ax-by+s,bx+ay+t) if sign=1,
        and (ax+by+s,bx-ay+t) if sign=-1
        '''
        if parent is None:
            raise ValueError("The parent must be provided")
        self._a=a
        self._b=b
        self._s=s
        self._t=t
        if not (sign==1 or sign==-1):
            raise ValueError("sign must be either 1 or -1.")
        self._sign=sign
        self._parent=parent
        MultiplicativeGroupElement.__init__(self,parent)

    def _mul_(self,s):
        r'''Compose two similarity reflections.'''
        C = self.__class__
        if self._sign==1:
            return C(self._parent,
                self._a*s._a-self._b*s._b, 
                self._b*s._a+self._a*s._b,
                self._a*s._s-self._b*s._t+self._s,
                self._b*s._s+self._a*s._t+self._t,s._sign)
        else:
            return C(self._parent,
                self._a*s._a+self._b*s._b, 
                self._b*s._a-self._a*s._b,
                self._a*s._s+self._b*s._t+self._s,
                self._b*s._s-self._a*s._t+self._t,-1*s._sign)

    def __invert__(self):
        r'''Invert a similarity.'''
        if self._sign==1:
            det=self._a*self._a+self._b*self._b
            a=self._a/det
            b=-self._b/det
            C = self.__class__
            return C(self._parent,
                a,
                b,
                -a*self._s+b*self._t,
                -b*self._s-a*self._t,1)
        else:
            det=-(self._a*self._a+self._b*self._b)
            a=-self._a/det
            b=-self._b/det
            C = self.__class__
            return C(self._parent,
                a,
                b,
                -a*self._s-b*self._t,
                -b*self._s+a*self._t,-1)

    def _div_(self,s):
        return self._mul_(s.__invert__())

    def __hash__(self):
        return 73*hash(self._a)-19*hash(self._b)+13*hash(self._s)+53*hash(self._t)+67*hash(self._sign)

    def __call__(self,w):
        r'''Return m*w+v.'''
        if self._sign==1:
            return vector([self._a*w[0]-self._b*w[1]+self._s, self._b*w[0]+self._a*w[1]+self._t])
        else:
            return vector([self._a*w[0]+self._b*w[1]+self._s, self._b*w[0]-self._a*w[1]+self._t])

    def _repr_(self):
        return "ReflectiveSimilarity (x,y) mapsto ("+str(self._a)+"*x-"+\
            str(self._sign*self._b)+"*y+"+str(self._s)+", "+\
            str(self._b)+"*x+"+str(self._sign*self._a)+"*y+"+str(self._t)+")"

    def _cmp_(self, other):
        x=cmp(self._a,other._a)
        if x!=0:
            return x
        x=cmp(self._b,other._b)
        if x!=0:
            return x
        x=cmp(self._s,other._s)
        if x!=0:
            return x
        x=cmp(self._t,other._t)
        if x!=0:
            return x
        return cmp(self._sign,other._sign)

    __cmp__=_cmp_

    def matrix(self):
        return matrix(self._parent._f,[
            [self._a, -self._sign*self._b, self._s],
            [self._b,  self._sign*self._a, self._t],
            [self._parent._f.zero(), self._parent._f.zero(), self._parent._f.one()]])

    def a(self):
        return self._a

    def b(self):
        return self._b

    def s(self):
        return self._s

    def t(self):
        return self._t

    def sign(self):
        return self._sign

    def derivative(self):
        r"""Return the 2x2 matrix corresponding to the derivative of the similarity of the plane."""
        return matrix(self._parent._f,[
            [self._a, -self._sign*self._b],
            [self._b,  self._sign*self._a]])


class SimilarityReflectionGroup(UniqueRepresentation,Group):
    r'''Group representing all possibly orientation reversing similarities in the plane.
    This is the group generated by rotations, translations and dilations.
    '''

    Element = SimilarityReflection

    def _element_constructor_(self, *args, **kwds):
        if len(args)!=1:
            return self.element_class(self, *args, **kwds)
        x = args[0]
        p=parent(x)
        if self._f.has_coerce_map_from(p):
            return self.element_class( self,self._f(x), self._f.zero(), self._f.zero(), self._f.zero())
        if isinstance(p, SimilarityReflectionGroup):
            return self.element_class(self, x.a(), x.b(), x.s(), x.t(), x.new())
        if isinstance(p, SimilarityGroup):
            return self.element_class(self, x.a(), x.b(), x.s(), x.t(), 1)
        if isinstance(p, TranslationGroup):
            return self.element_class( self,self._f.one(), self._f.zero(), x.s(), x.t(), 1)
        return self.element_class(self, x, **kwds)

    def _coerce_map_from_(self, S):
        if self._f.has_coerce_map_from(S):
            return True
        if isinstance(S, SimilarityReflectionGroup):
            return self._f.has_coerce_map_from(S._f)
        if isinstance(S, SimilarityGroup):
            return self._f.has_coerce_map_from(S._f)
        if isinstance(S, TranslationGroup):
            return self._f.has_coerce_map_from(S.base_field())
           
    def __init__(self, base_field):
        self._f=base_field
        # The vector space of vectors 
        self._vs = VectorSpace(self._f,2)
        Group.__init__(self, category=Groups().Infinite())

    def _repr_(self):
        return "SimilarityReflectionGroup over field "+str(self._f)

    def one(self):
        return self.element_class(self,self._f.one(),self._f.zero(),self._f.zero(),self._f.zero(),1)

    def an_element(self):
        return self.element_class(self,self._f(ZZ_3),self._f(ZZ_4),self._f(ZZ_2),self._f(-ZZ_1),-1)

    def is_abelian(self):
        return False

    def gens(self):
        pairs=[
            (self._f.one(),self._f.zero()),
            (self._f(ZZ_2),self._f.zero()),
            (self._f.zero(),self._f(ZZ_2)),
            (self._f(ZZ_3),self._f(ZZ_4))]
        l=[]
        for p in pairs:
            for v in self._vs.gens():
                l.append(self.element_class(self,p[0],p[1],v[0],v[1],1))
                l.append(self.element_class(self,p[0],p[1],v[0],v[1],-1))
        return l

    def base_field(self):
        return self._f


