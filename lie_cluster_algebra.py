from sage.algebras.cluster_algebra import ClusterAlgebra

class CoxeterClusterAlgebra(ClusterAlgebra):

    @staticmethod
    def __classcall__(self, ct, **kwargs):
        kwargs['cartan_type'] = CartanType(ct)
        n = kwargs['cartan_type'].rank()
        kwargs.setdefault('coxeter_element', range(n))
        kwargs['coxeter_element'] = tuple(kwargs['coxeter_element'])
        c = { i:j for (j,i) in enumerate(kwargs['coxeter_element']) }
        B = kwargs['cartan_type'].cartan_matrix() - 2
        for i in range(n):        
            for j in range(n):    
                if c[i] < c[j]:         
                    B[i,j] = -B[i,j] 
        B = block_matrix([[B],[identity_matrix(n)],[identity_matrix(n)]])
        return super(CoxeterClusterAlgebra, self).__classcall__(self, B, **kwargs)

    def __init__(self, B, **kwargs):
        ClusterAlgebra.__init__(self, B, **kwargs)
        self._cartan_type = kwargs['cartan_type']
        self._coxeter_element = kwargs['coxeter_element']
        self._subs = { self.initial_cluster_variable(i):var('h%d'%i) for i in range(B.ncols()) }
       
    @cached_method
    def cluster_variable(self, g_vector):
        x = super(CoxeterClusterAlgebra, self).cluster_variable._instance_call(g_vector)
        return x.lift().subs(self._subs)

class OldFoo(UniqueRepresentation):

    @staticmethod
    def __classcall__(self, data, **kwargs):
        hashable_data = tuple(data)
        kwargs['some_default_option'] = 'bar'
        return super(OldFoo, self).__classcall__(self, hashable_data, **kwargs)

    def __init__(self, data, **kwargs):
        self._data = data
        self._some_default_option = kwargs['some_default_option']

    def __repr__(self):
        return 'Data = %s, option = %s'%(self._data, self._some_default_option)

class NewFoo(OldFoo):

    @staticmethod
    def __classcall__(self, data, **kwargs):
        new_data = data+data
        return super(NewFoo, self).__classcall__(self, new_data, **kwargs)

    def __init__(self, data, **kwargs):
        OldFoo.__init__(self, data, **kwargs)
        self._some_default_option = 'foobar'

