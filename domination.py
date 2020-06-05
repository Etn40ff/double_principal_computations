class DominanceOrder(SageObject):
    r"""
    A class to compute the dominace order as defined in `arXiv:1902.09507`
    """

    def __init__(self, B):
        self.B = copy(B)
        self.A = ClusterAlgebra(B,principal_coefficients=True)

    def paths_up_to_length(self, k):
        paths = [ s.path_from_initial_seed() for s in self.A.seeds(mutating_F=False, depth=k) ]
        prefixes = [ p[:-1] for p in paths ]
        return [ p for p in paths if p not in prefixes ]

    @cached_method(key=lambda self, g, depth: (tuple(g), depth))
    def dominated_polytope(self, g, depth):
        paths = self.paths_up_to_length(depth)
        p = paths.pop()
        polytope = cut_along_sequence(g,self.B,p)
        while paths:
            p = paths.pop()
            polytope = polytope.intersection(cut_along_sequence(g,self.B,p))
        return polytope

    def dominated_g_vectors(self, g, depth):
        polytope = self.dominated_polytope(g, depth)
        pts = polytope.integral_points()
        lattice = self.B.transpose().image()
        g = vector(g)
        return [ p for p in pts if p-g in lattice ]

    def show_domination(self, g, depth):
        g = vector(g)
        pts = self.dominated_g_vectors(g, depth)
        polytope = self.dominated_polytope(g, depth)
        return polytope.plot(zorder=-1) + list_plot(pts, color='red', size=50) 

def cut_along_sequence(g, B, seq):
    current_polytope = Polyhedron(rays=B.columns(),base_ring=QQ).translation(g)
    if seq == []:
        return current_polytope

    k = seq.pop()
    n = B.ncols()
    I = identity_matrix(n)

    Hp = Polyhedron(ieqs=[(0,)*(k+1)+(1,)+(0,)*(n-k-1)])
    Mp = copy(I)
    Mp[:,k] = matrix(n,1, lambda i,_: max(B[i,k],0) if i != k else -1)

    Hm = Polyhedron(ieqs=[(0,)*(k+1)+(-1,)+(0,)*(n-k-1)])
    Mm = copy(I)
    Mm[:,k] = matrix(n,1, lambda i,_: max(-B[i,k],0) if i != k else -1)
    
    new_g = (Mp if g in Hp else Mm) * vector(g)
    new_B = copy(B)
    new_B.mutate(k)

    new_polytope = cut_along_sequence(new_g, new_B, seq)
    
    new_polytope_p = Mm*(new_polytope.intersection(Hp))
    new_polytope_m = Mp*(new_polytope.intersection(Hm))
    
    new_polytope = new_polytope_p.convex_hull(new_polytope_m)

    return new_polytope.intersection(current_polytope)