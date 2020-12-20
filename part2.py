import numpy as np
import scipy.sparse as sp
import time
import matplotlib.pyplot as plt

def graphGen(n,a):
    K = np.round(n**(1/(a-1)))
    k = np.arange(K)+1
    p = k**(-a)
    p = p/np.sum(p)
    cmf = np.cumsum(p)
    degrees = np.zeros(n,dtype=int)
    for i in range(n):
        degrees[i] = 1+int(len(np.argwhere(np.random.rand()>cmf).flatten()))
    num_stubs = np.sum(degrees)
    if num_stubs%2==1:
        j = np.random.randint(0,n)
        degrees[j]=degrees[j]+1
        num_stubs = num_stubs+1
    num_edges = int(num_stubs/2)
    bag = np.array([])
    for i in range(n):
        aux = i*np.ones(degrees[i],dtype=int)
        bag = np.append(bag,aux)
    bag = np.random.permutation(bag)
    bag1 = bag[:num_edges]
    bag2 = bag[num_edges:]
    edges = np.zeros((n,n))
    for i in range(num_edges):
        a = int(bag1[i])
        b = int(bag2[i])
        if a!=b:
            edges[a,b]=1
            edges[b,a]=1
    return edges

def transmissionMatrix(mat,P):
    m = mat.copy()
    inds = m.nonzero()
    l = len(inds[0])
    for i in range(l):
        if np.random.rand()>P:m[inds[0][i],inds[1][i]]=0
    return m

class Node:
    def __init__(self, index, nodes):
        self.num = index
        self.label = -1
        self.dist = np.inf
        self.parent = -1
        self.td = np.inf
        self.tf = np.inf
        self.neighbors = nodes
        
class Graph:
    def __init__(self, mat):
        self.size = mat.shape[0]
        self.nodes = [Node(i,mat[i].nonzero()[0]) for i in range(self.size)]
        self.time = 0
        
    def BFS(self, s):
        for node in self.nodes:
            node.parent = -1
            node.dist = np.inf
            node.label = -1
        
        self.nodes[s].label+=1
        self.nodes[s].dist=0
                       
        Q = [self.nodes[s]]
        
        while Q:
            u = Q[0]
            ns = u.neighbors
            for n in ns:
                v = self.nodes[n]
                if v.label == -1:
                    v.label+=1
                    v.dist = u.dist+1
                    v.parent = u.num
                    Q.append(v)
            u.label+=1
            Q.pop(0)
            
        comp = np.array([ j for j in range(self.size) if self.nodes[j].label==1 ])
        dists = np.array([ self.nodes[i].dist for i in comp ])
        parents = np.array([ self.nodes[i].parent for i in comp ])        
        return (comp,dists,parents)
    
    def matrix(self):
        m = np.zeros((self.size,self.size))
        for i in range(self.size):
            m[i] = np.array([ a in self.nodes[i].neighbors for a in range(self.size)]).astype(int)
        return m
    
    def DFS(self):
        comps = []
        incomplete = True
        unvisited = np.array([ i for i in range(self.size) if self.nodes[i].label==-1 ])
        visited = []
        while incomplete:
            u = self.nodes[unvisited[0]]
            self.DFS_Visit(u)
            unvisited = np.array([ i for i in range(self.size) if self.nodes[i].label==-1 ])
            visited = np.setdiff1d(np.arange(self.size),unvisited)
            for c in comps:
                visited = np.setdiff1d(visited,c)
            comps.append(visited)
            incomplete = len(unvisited)>0
        return comps
            
    
    def DFS_Visit(self, u):
        self.time+=1
        u.td=self.time
        u.label+=1
        for n in u.neighbors:
            v = self.nodes[n]
            if v.label==-1:
                v.parent=u.num
                self.DFS_Visit(v)
        u.label+=1
        self.time+=1
        u.tf=self.time
        
# do parts 3 a-b
t0 =time.time()
print(t0)
a = 2.2
n=10000
T = 0.4
runs = 20
giantS = np.zeros(runs)
epidemicS = np.zeros(runs)

for i in range(runs):
    mat = graphGen(n,a)    
    g = Graph(mat)    
    comps = g.DFS()
    sizes = np.array([ len(c) for c in comps ])
    giantS[i]=np.max(sizes)/n
    t1 = time.time()
    print(t1-t0)
    
    m2 = transmissionMatrix(mat,T)    
    g2 = Graph(m2)
    comps2 = g2.DFS()
    sizes2 = np.array([ len(c) for c in comps2 ])
    epidemicS[i] = np.max(sizes2)/n

print(time.time()-t1)
print('average giant size'+str(np.mean(giantS)))
print('average epidemic size'+str(np.mean(epidemicS)))

# do 3c

#Ts = np.linspace(0.025,0.2,8)
#runs = 5
#avEpiS = np.zeros(8)
#for i,t in enumerate(Ts):
#    aux = np.zeros(runs)
#    for j in range(runs):
#        mat = graphGen(n,a)
#        m2 = transmissionMatrix(mat,t)
#        g = Graph(m2)
#        comps = g.DFS()
#        sizes = np.array([ len(c) for c in comps2 ])
#        aux[j] = np.max(sizes)/n
#    avEpiS[i]=np.mean(aux)
#
#fig=plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(Ts,avEpiS)
#ax.set_xlabel('T')
#ax.set_ylabel('Average Fraction Infected')
#plt.show()


## do part 4

#t0=time.time()
#a = 2.2
#n=10000
#T = 0.4
#runs = 100
#starts = np.random.permutation(n)[:runs]
#mat = graphGen(n,a)
#m2 = transmissionMatrix(mat,T)
#g = Graph(m2)
#n_infs = []
#for i in range(runs):
#    (comp,dists,parents) = g.BFS(starts[i])
#    d = np.max(dists)
#    n_inf = np.zeros(d+1)
#    for j in range(d+1):
#        n_inf[j]= len(np.argwhere(dists==j).flatten())/n
#    n_infs.append(n_inf)
#
#fig=plt.figure()
#ax = fig.add_subplot(111)
#for s in n_infs:
#    ax.plot(s)
#ax.set_xlabel('Time')
#ax.set_ylabel('Number Infected')
#plt.show()
#print(time.time()-t0)    