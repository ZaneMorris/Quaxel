##############################################################################
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D as Line

class point:
    def __init__(self,pt):
        self.x, self.y = pt[0], pt[1]

class node:
    def __init__(self, center, el, z_el, lvl=0, tris = None):
        self.c = np.array(center)
        self.e0 = np.array([1,0,0])
        self.e1 = np.array([0,1,0])
        self.e2 = np.array([0,0,1])
        self.el, self.z_el = el, z_el
        self.tris = []
        self.lvl = lvl
        self.leaf_flag = False
    
    def intersect(self,tri):
        #returns true if there exists collision - no axis that seperates the two (no collision)
        #returns false if no collision exists - there is a seperating axes

        #translate triangle vertices toward origin, s.t. self is located at origin
        v0 = np.array(tri.v0.val) - self.c #element-wise subtraction
        v1 = np.array(tri.v1.val) - self.c #element-wise subtraction
        v2 = np.array(tri.v2.val) - self.c #element-wise subtraction
        tri_edges = [v1-v0, v2-v1, v0-v2]

        #Test Axes - potential seperating axes
        test_axes = [self.e0,self.e1,self.e2] #basis vectors of the axes the self is aligned with
        edge_axes = [np.cross(e,f)/np.linalg.norm(np.cross(e,f),ord=2) for f in tri_edges for e in test_axes] #mutally orthog. vects for bases and triangle edges
        test_axes.extend(edge_axes) #entire list of test axes, short of facet normal
        
        #prep
        x = self.el/2
        z = self.z_el/2

        #First Test Set
        for a in test_axes:
            r = x*abs(a[0]) + x*abs(a[1]) + z*abs(a[2]) #effective radius of projection for self
            p0, p1, p2 = np.dot(a,v0), np.dot(a,v1), np.dot(a,v2)
            if min(p0,p1,p2) > r or max(p0,p1,p2) < -r:
                #print('no collision')
                return False
        
        #Second Test Set
        vertices = [[x,x,z],[x,-x,z],[-x,-x,z],[-x,x,z],[x,x,-z],[x,-x,-z],[-x,-x,-z],[-x,x,-z]] #raw, origin-centered AABB vertices
        vertices = [(np.array(v)+self.c)-tri.cog for v in vertices] #translating vector root to facet origin/cog
        val0 = np.dot(tri.n,vertices[0])

        for i in range(1,len(vertices)):
            val = np.dot(tri.n,vertices[i])
            if val*val0 <= 0:
                break
            elif i==7:
                #print('no collision')
                return False #all vertices are on 1 side of plane, therefore a seperating axis
        
        #print('collision')
        return True #collision exists

    def get_intersections(self, parent_tris):
        for tri in parent_tris:
            if self.intersect(tri):
                self.tris.append(tri)

    def subdivide(self, th_nd, th_br):
        if len(self.tris) <= th_nd or self.lvl >= th_br:
            self.leaf_flag = True
            return
        else:
            self.leaf_flag = False #here for add_point() functionality
        
        self.child1 = node([self.c[0]+self.el/4, self.c[1]+self.el/4, self.c[2]], self.el/2, self.z_el, self.lvl+1)
        self.child1.get_intersections(self.tris)
        self.child1.subdivide(th_nd, th_br)

        self.child2 = node([self.c[0]-self.el/4, self.c[1]+self.el/4, self.c[2]], self.el/2, self.z_el, self.lvl+1)
        self.child2.get_intersections(self.tris)
        self.child2.subdivide(th_nd, th_br)

        self.child3 = node([self.c[0]-self.el/4, self.c[1]-self.el/4, self.c[2]], self.el/2, self.z_el, self.lvl+1)
        self.child3.get_intersections(self.tris)
        self.child3.subdivide(th_nd, th_br)

        self.child4 = node([self.c[0]+self.el/4, self.c[1]-self.el/4, self.c[2]], self.el/2, self.z_el, self.lvl+1)
        self.child4.get_intersections(self.tris)
        self.child4.subdivide(th_nd, th_br)

    def draw_node(self,lines=None):
        if lines is None:
            lines = []
        
        x_1, x_2 = self.c[0] - self.el/2, self.c[0] + self.el/2
        y_1, y_2 = self.c[1] - self.el/2, self.c[1] + self.el/2
        lines.extend([Line([x_1,x_1],[y_1,y_2]), Line([x_1,x_2],[y_2,y_2]), Line([x_2,x_2],[y_2,y_1]), Line([x_2,x_1],[y_1,y_1])])
        
        if not self.leaf_flag:
            lines = self.child1.draw_node(lines)
            lines = self.child2.draw_node(lines)
            lines = self.child3.draw_node(lines)
            lines = self.child4.draw_node(lines)
        
        return lines

class qtree:
    
    def __init__(self, facets, th_nd=2, th_br = 6, boundary = []):
        """boundary: [x_min, x_max, y_min, y_max, z_min, z_max]"""
        
        #Make boundary values square around center, with a proportionally tailored offset
        x_r, y_r = abs(boundary[1]-boundary[0]), abs(boundary[3]-boundary[2])
        max_r = max(x_r,y_r)
        ofs = max_r/100
        boundary[0] = boundary[0] - (max_r-x_r)/2 - ofs
        boundary[1] = boundary[1] + (max_r-x_r)/2 + ofs
        boundary[2] = boundary[2] - (max_r-y_r)/2 - ofs
        boundary[3] = boundary[3] + (max_r-y_r)/2 + ofs

        #Attribute Definition
        self.boundary = boundary
        self.th_nd, self.th_br = th_nd, th_br
        self.center = [(boundary[0] + boundary[1]) / 2, (boundary[2] + boundary[3]) / 2, (boundary[4] + boundary[5]) / 2]
        self.root = node(self.center, abs(boundary[1]-boundary[0]), abs(boundary[5]-boundary[4]))
        
        #Initial Intersections (Should be all facets)
        self.root.get_intersections(facets)

    def subdivide(self):
        self.root.subdivide(self.th_nd, self.th_br)

    def find_leaf(self, point):
        node = self.root
        while not node.leaf_flag:
            if point[1] > node.c[1]:
                if point[0] > node.c[0]:
                    node = node.child1
                else:
                    node = node.child2
            else:
                if point[0] > node.c[0]:
                    node = node.child4
                else:
                    node = node.child3
        return node.tris

    def voxel_mesh(self,n):
        #given bounds and number of voxels 
        el = min(abs(self.boundary[1]-self.boundary[0]), abs(self.boundary[3]-self.boundary[2]), abs(self.boundary[5]-self.boundary[4]))/n
        x_n = int(np.ceil(abs(self.boundary[1]-self.boundary[0])/el))
        y_n = int(np.ceil(abs(self.boundary[3]-self.boundary[2])/el))
        z_n = int(np.ceil(abs(self.boundary[5]-self.boundary[4])/el))

        self.vox_arr = np.zeros((z_n,x_n,y_n)) #by default, all values are 0
        return el

    def voxelize(self, el):
        x = self.boundary[0]+el/2
        for i in range(self.vox_arr.shape[1]): #in this case, index 1 is the x dim
            y = self.boundary[2]+el/2
            for j in range(self.vox_arr.shape[2]): #in this case, index 2 is the y dim
                tris = self.find_leaf([x,y]) #list of intersecting facets
                r1, r2 = (el/16)*np.random.uniform(low=-1.0,high=1.0), (el/8)*np.random.uniform(low=-1.0,high=1.0)
                #r1 and r2 supposed to hypothetically alleviate the chances of aligning ray with edge shared by 2 facets
                origin = [x+r1,y+r2,self.boundary[4]]
                z_col = ray_intersection(self.vox_arr, origin, el, tris)
                
                #update vox_arr
                self.vox_arr[:,i,j] = z_col
                y += el
            x += el
        return

    def FEM(self):
        
        return

    def draw_quad(self,**kwargs):
        lines = self.root.draw_node()
        fig, ax = plt.subplots(**kwargs)
        ax.set_xlim(left=self.boundary[0],right=self.boundary[1])
        ax.set_ylim(bottom=self.boundary[3],top=self.boundary[2])
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')
        for line in lines:
            ax.add_line(line)
        plt.show()

    def draw_temp(self, el, **kwargs):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        x = self.boundary[0]+(el/2)

        for i in range(self.vox_arr.shape[1]):
            y = self.boundary[2]+(el/2)
            for j in range(self.vox_arr.shape[2]):
                z = self.boundary[4]+(el/2)
                for k in range(self.vox_arr.shape[0]):
                    if self.vox_arr[k,i,j]:
                        pl = np.array([x,y,z])+(el/2)
                        ne = np.array([x,y,z])-(el/2)
                        
                        X,Y = np.meshgrid([ne[0],pl[0]],[ne[1],pl[1]])
                        ax.plot_surface(X,Y,ne[2]*np.ones((2,2)), **kwargs) #low face
                        ax.plot_surface(X,Y,pl[2]*np.ones((2,2)), **kwargs) #low face
                        
                        X,Z = np.meshgrid([ne[0],pl[0]],[ne[2],pl[2]])
                        ax.plot_surface(X,ne[1]*np.ones((2,2)),Z, **kwargs) #front face
                        ax.plot_surface(X,pl[1]*np.ones((2,2)),Z, **kwargs) #front face

                        Y,Z = np.meshgrid([ne[1],pl[1]],[ne[2],pl[2]])
                        ax.plot_surface(ne[0]*np.ones((2,2)),Y,Z, **kwargs) #right face
                        ax.plot_surface(pl[0]*np.ones((2,2)),Y,Z, **kwargs) #left face

                    z += el
                y += el
            x += el
        ax.set_xlabel('x [mm]')
        ax.set_ylabel('y [mm]')
        ax.set_zlabel('z [mm]')
        plt.show()
        return

    def draw2(self):
        #axis length ratios
        x_rn, y_rn, z_rn = abs(self.boundary[1]-self.boundary[0]), abs(self.boundary[3]-self.boundary[2]), abs(self.boundary[5]-self.boundary[4])
        max_rn = max(x_rn,y_rn,z_rn)
        xr, yr, zr = x_rn/max_rn, y_rn/max_rn, z_rn/max_rn

        #z,x,y = np.indices((self.vox_arr.shape[0],self.vox_arr.shape[1],self.vox_arr.shape[2]))
        ax = plt.figure().add_subplot(projection='3d')
        ax.voxels(self.vox_arr, facecolors='b', edgecolor='k')
        ax.set_box_aspect([zr,xr,yr])
        ax.set_xlabel('z [mm]')
        ax.set_ylabel('x [mm]')
        ax.set_zlabel('y [mm]')
        plt.show()
        return

################## FUNCTIONS ##################

def ray_intersection(vox_arr,origin,el,tris):
    z_col = np.zeros((vox_arr.shape[0])) #vox_arr.shape[0] is len of array in z axis
    
    #default: bottom voxel not filled (vox_arr)
    z_inds = [] #z locations along ray that intersect with facet
    r = np.array([0,0,1]) #ray direction vector
    
    for tri in tris:
        #distance along facet normal from origin [0,0,0] to plane
        d = -1*(tri.n[0]*tri.v0.val[0] + tri.n[1]*tri.v0.val[1] + tri.n[2]*tri.v0.val[2])
        
        # checking perpindicularity of facet normal and ray
        denom = np.dot(tri.n,r)
        if denom > 1e-10 or denom < -1e-10:
            #find distance along r to intersection with facet plane
            t = -1*(np.dot(tri.n,origin)+d) / denom
            p = origin + t*r
        else:
            #ray direction and facet normal are perpindicular, no intersection
            continue #return to top of for loop
        
        #Inside-Outside Test
        #facet vertices are listed in ccw order when looking along normal from outside
        v10 = tri.v1.val - tri.v0.val
        v21 = tri.v2.val - tri.v1.val
        v02 = tri.v0.val - tri.v2.val
        if np.dot(tri.n,np.cross(v10,p-tri.v0.val)) > 0 and \
           np.dot(tri.n,np.cross(v21,p-tri.v1.val)) > 0 and \
           np.dot(tri.n,np.cross(v02,p-tri.v2.val)) > 0:
            #lies on plane
            z_inds.append(round(p[2]/el)) #this value should be the cap for next
        else:
            #lies outside the confines of the facet, therefore no intersection
            continue #return to top of for loop
    
    val = 0
    ind0 = 0
    z_inds.sort()
    for i in range(len(z_inds)):
        if ind0 == z_inds[i]: #end condition and "cantilever" condition
            z_col[ind0-1] = val
        else:
            z_col[ind0:z_inds[i]] = val
        ind0 = z_inds[i]
        val = 0 if val else 1
    return z_col
