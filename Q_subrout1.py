import numpy as np

##############################################################################
# Contains classes: Vertex, Facet
# Contains functions: read_stl, update_bounds, utls

def utl0(vert):
    return vert[0]
def utl1(vert):
    return vert[1]
def utl2(vert):
    return vert[2]

class vertex:
    def __init__(self,v,flag):
        self.val, self.flag = np.array(v), flag

class facet:
    def __init__(self,v0,v1,v2,n):
        self.v0, self.v1, self.v2 = vertex(v0,0), vertex(v1,1), vertex(v2,2)
        self.n = n
        x = (self.v0.val[0]+self.v1.val[0]+self.v2.val[0])/3
        y = (self.v0.val[1]+self.v1.val[1]+self.v2.val[1])/3
        z = (self.v0.val[2]+self.v1.val[2]+self.v2.val[2])/3
        self.cog = [x,y,z]
        #self.v0flag, self.v1flag, self.v2flag = 0, 1, 2

def read_stl(path,filename):
    """
    Returns a list of facet objects, with v0, v1, v2 attributes, themselves vertex objects, with
    a .val attribute, which is a list of [x,y,z] coordinates
    """
    f = open(path+'\\'+filename,'r')
    lines = f.readlines()
    f.close()

    facets = []
    bounds = [1e6, -1e6, 1e6, -1e6, 1e6, -1e6] #[x_min, x_max, y_min, y_max, z_min, z_max]
    #ind = 3
    ind = -4
    while lines[ind+5] != 'endsolid':
        ind += 7
        n = lines[ind-2][16:-1].split(' ')
        n = np.array([float(n[i]) for i in range(3)])
        v0 = lines[ind][16:-1].split(' ')
        v0 = [float(v0[i]) for i in range(3)]
        v1 = lines[ind+1][16:-1].split(' ')
        v1 = [float(v1[i]) for i in range(3)]
        v2 = lines[ind+2][16:-1].split(' ')
        v2 = [float(v2[i]) for i in range(3)]
        facets.append(facet(v0,v1,v2,n))
        #ind += 7

        # update bounds
        bounds = update_bounds(bounds, [v0,v1,v2])


    del lines
    return facets, bounds

def update_bounds(bounds, vertices):
    bounds[0] = min(min(vertices, key=utl0)[0], bounds[0])
    bounds[1] = max(max(vertices, key=utl0)[0], bounds[1])
    bounds[2] = min(min(vertices, key=utl1)[1], bounds[2])
    bounds[3] = max(max(vertices, key=utl1)[1], bounds[3])
    bounds[4] = min(min(vertices, key=utl2)[2], bounds[4])
    bounds[5] = max(max(vertices, key=utl2)[2], bounds[5])
    return bounds