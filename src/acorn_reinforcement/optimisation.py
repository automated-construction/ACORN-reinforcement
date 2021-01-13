from math import gcd, ceil

import itertools as it

import numpy as np
from scipy import sparse
import cvxpy as cvx

import matplotlib.pyplot as plt


__all__ = ['lo', 'lo2']


def lo2(nodes, lines, joint_cost=0.0):

    #Calculate equilibrium matrix B
    def calcB(Nd, Cn, dof):
        m, n1, n2 = len(Cn), Cn[:,0].astype(int), Cn[:,1].astype(int)
        l, X, Y, Z = Cn[:,2], Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = dof[n1*3], dof[n1*3+1], dof[n1*3+2], dof[n2*3], dof[n2*3+1], dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l * d2, X/l * d3, Y/l * d4, Z/l *d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        return sparse.coo_matrix((s, (r, c)), shape = (len(Nd)*3, m))
    #Solve linear programming problem
    def solveLP(Nd, Cn, f, dof, jc):
        l = [col[2] + jc for col in Cn]
        B = calcB(Nd, Cn, dof)
        a = cvx.Variable(len(Cn))
        obj = cvx.Minimize(np.transpose(l) @ a)
        q, eqn, cons = [],  [], [a>=0]
        for k, fk in enumerate(f):
            q.append(cvx.Variable(len(Cn)))
            eqn.append(B * q[k] == fk * dof)
            cons.extend([eqn[k], q[k] >= - cvx.multiply(Cn[:,5], a), q[k] <= cvx.multiply(Cn[:,4], a)])
        prob = cvx.Problem(obj, cons)
        vol = prob.solve()
        q = [np.array(qi.value).flatten() for qi in q]
        a = np.array(a.value).flatten()
        u = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]
        return vol, a, q, u
    #Check dual violation
    def stopViolation(Nd, PML, dof, u, jc):
        lst = np.where(PML[:,3]==False)[0]
        Cn = PML[lst]
        l = Cn[:,2] + jc
        B = calcB(Nd, Cn, dof).tocsc()
        y = np.zeros(len(Cn))
        for uk in u:
            yk = np.multiply(B.transpose().dot(uk) / l, np.array([Cn[:,4], -Cn[:,5]]))
            y += np.amax(yk, axis=0)
        vioCn = np.where(y>1.0001)[0]
        vioSort = np.flipud(np.argsort(y[vioCn]))
        num = ceil(min(len(vioSort), 0.05*max( [len(Cn)*0.05, len(vioSort)])))
        for i in range(num): 
            PML[lst[vioCn[vioSort[i]]]][3] = True
        return num == 0

    supports = {}
    loads = [{} for _ in range((len(nodes[0]) - 2))]
    for key, data in enumerate(nodes):
        supports[key] = data[1]
        for i, d in enumerate(data[2:]):
            loads[i][key] = d
    nodes = [data[0] for data in nodes]

    # points
    Nd = np.array(nodes)
    dof, f, PML = np.ones((len(Nd),3)), [], []
    # supports
    for i, supp in supports.items():
        dof[i,:] = supp
    # loads
    for loads_k in loads:
        for i, nd in enumerate(Nd):
            f += loads_k[i] if i in loads_k else [0, 0, 0]
    # ground structure
    for i, j, inc, st, sc in lines:
        l = np.sqrt(sum([(Nd[i][xyz] - Nd[j][xyz])**2 for xyz in range(3)]))
        PML.append( [i, j, l, inc, st, sc] )
    PML, dof = np.array(PML), np.array(dof).flatten()
    f = [f[i:i+len(Nd)*3] for i in range(0, len(f), len(Nd)*3)]
    # start adding members
    for itr in range(1, 100):
        Cn = PML[PML[:,3] == True]
        vol, a, q, u = solveLP(Nd, Cn, f, dof, joint_cost)
        if stopViolation(Nd, PML, dof, u, joint_cost): break
    return Nd, Cn, a, q

def lo(nodes, lines, supports, loads, joint_cost=0.0):
    """Layout optimisation based on He et al., 2019.

    Parameters
    ----------
    nodes : list
        Node data as list of XYZ coordinates.
    lines : list
        Line data as list of tuples (start node index, end node index, include in starting set, tensile strength, compressive strength).
        Strengths are assumed to be in MPa.
        The line is included in the initial set if the boolean is True.
    supports : dict
        Support data as a dictionary of node indices to support conditions.
        The support condition is a list of 3 binaries for each X, Y, Z direction, with 0 meaning fixed and 1 free.
        Nodes with no specifications are free ([1, 1, 1]).
    loads : list
        Load data as a list of dictionary of node indices to load vectors.
        Each list item represents a load case.
        The load vector is a list of 3 floats for each X, Y, Z direction.
        Nodes with no specifications are not loaded ([0.0, 0.0, 0.0]).
        Loads are assumed to be in kN.
    joint_cost : float
        Cost of the joints. Added as extra length on members.
        Defaults is 0.0.

    Returns
    -------
    Nd: Numpy array of node vertices.
    Cn: Numpy array of used lines with i, j, l, inc, st, sc.
    a: Numpy array of cross-sections.
    q: Numpy array of forces.

    Reference
    ---------
    L. He, M. Gilbert, X. Song, "A Python script for adaptive layout optimization of trusses", Struct. Multidisc. Optim., 2019. 

    """

    #Calculate equilibrium matrix B
    def calcB(Nd, Cn, dof):
        m, n1, n2 = len(Cn), Cn[:,0].astype(int), Cn[:,1].astype(int)
        l, X, Y, Z = Cn[:,2], Nd[n2,0]-Nd[n1,0], Nd[n2,1]-Nd[n1,1], Nd[n2,2]-Nd[n1,2]
        d0, d1, d2, d3, d4, d5 = dof[n1*3], dof[n1*3+1], dof[n1*3+2], dof[n2*3], dof[n2*3+1], dof[n2*3+2]
        s = np.concatenate((-X/l * d0, -Y/l * d1, -Z/l * d2, X/l * d3, Y/l * d4, Z/l *d5))
        r = np.concatenate((n1*3, n1*3+1, n1*3+2, n2*3, n2*3+1, n2*3+2))
        c = np.concatenate((np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m), np.arange(m)))
        return sparse.coo_matrix((s, (r, c)), shape = (len(Nd)*3, m))
    #Solve linear programming problem
    def solveLP(Nd, Cn, f, dof, jc):
        l = [col[2] + jc for col in Cn]
        B = calcB(Nd, Cn, dof)
        a = cvx.Variable(len(Cn))
        obj = cvx.Minimize(np.transpose(l) @ a)
        q, eqn, cons = [],  [], [a>=0]
        for k, fk in enumerate(f):
            q.append(cvx.Variable(len(Cn)))
            #print(B)
            #print(q[k])
            #print(f[k])
            #print(dof)
            eqn.append(B * q[k] == fk * dof)
            cons.extend([eqn[k], q[k] >= - cvx.multiply(Cn[:,5], a), q[k] <= cvx.multiply(Cn[:,4], a)])
        prob = cvx.Problem(obj, cons)
        vol = prob.solve()
        q = [np.array(qi.value).flatten() for qi in q]
        a = np.array(a.value).flatten()
        u = [-np.array(eqnk.dual_value).flatten() for eqnk in eqn]
        return vol, a, q, u
    #Check dual violation
    def stopViolation(Nd, PML, dof, u, jc):
        lst = np.where(PML[:,3]==False)[0]
        Cn = PML[lst]
        l = Cn[:,2] + jc
        B = calcB(Nd, Cn, dof).tocsc()
        y = np.zeros(len(Cn))
        for uk in u:
            yk = np.multiply(B.transpose().dot(uk) / l, np.array([Cn[:,4], -Cn[:,5]]))
            y += np.amax(yk, axis=0)
        vioCn = np.where(y>1.0001)[0]
        vioSort = np.flipud(np.argsort(y[vioCn]))
        num = ceil(min(len(vioSort), 0.05*max( [len(Cn)*0.05, len(vioSort)])))
        for i in range(num): 
            PML[lst[vioCn[vioSort[i]]]][3] = True
        return num == 0

    # points
    Nd = np.array(nodes)
    dof, f, PML = np.ones((len(Nd),3)), [], []
    # supports
    for i, supp in supports.items():
        dof[i,:] = supp
    # loads
    for loads_k in loads:
        for i, nd in enumerate(Nd):
            f += loads_k[i] if i in loads_k else [0, 0, 0]
    # ground structure
    for i, j, inc, st, sc in lines:
        l = np.sqrt(sum([(Nd[i][xyz] - Nd[j][xyz])**2 for xyz in range(3)]))
        PML.append( [i, j, l, inc, st, sc] )
    PML, dof = np.array(PML), np.array(dof).flatten()
    f = [f[i:i+len(Nd)*3] for i in range(0, len(f), len(Nd)*3)]
    # start adding members
    for itr in range(1, 100):
        Cn = PML[PML[:,3] == True]
        vol, a, q, u = solveLP(Nd, Cn, f, dof, joint_cost)
        if stopViolation(Nd, PML, dof, u, joint_cost): break
    return Nd, Cn, a, q


def lo_plot(Nd, Cn, a, q, threshold, str, update=True, plane='xy'):
    plt.ion() if update else plt.ioff()
    plt.clf(); plt.axis('off'); plt.axis('equal');  plt.draw()
    plt.title(str)
    tk = 5 / max(a)
    for i in [i for i in range(len(a)) if a[i] >= threshold]:
        if all([q[lc][i]>=0 for lc in range(len(q))]): c = 'r'
        elif all([q[lc][i]<=0 for lc in range(len(q))]): c = 'b'
        else: c = 'tab:gray'
        pos = Nd[Cn[i, [0, 1]].astype(int), :]
        if plane == 'xz':
            u, v = 0, 2
        elif plane == 'yz':
            u, v = 1, 2
        else:
            u,v = 0, 1
        plt.plot(pos[:, u], pos[:, v], c, linewidth = a[i] * tk) # choose plan for plot
    plt.pause(0.01) if update else plt.show()


if __name__ =='__main__':
    
    import itertools as it

    nodes = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
    lines = [(i, j, 1, 1, 1) for i, j in it.combinations(range(len(nodes)), 2)]
    supports = {0: [0, 0, 0], 1: [1, 0, 0]}
    loads = [{2: [10.0, 0.0, 0.0]}]
    #loads = [{2: [1.0, 0.0, 0.0]}, {3: [1.0, 1.0, 0.0]}]

    Nd, Cn, a, q = lo(nodes, lines, supports, loads)
    lo_plot(Nd, Cn, a, q, max(a) * 1e-3, str='Optimal layout', update=False, plane='xy')
    print(a)

    # nodes = [
    # ((0.0, 0.0, 0.0), [0, 0, 0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    # ((1.0, 0.0, 0.0), [1, 0, 0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    # ((1.0, 1.0, 0.0), [1, 1, 1], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
    # ((0.0, 1.0, 0.0), [1, 1, 1], [0.0, 0.0, 0.0], [1.0, 1.0, 0.0]),
    # ]
    # lines = [(i, j, 1, 1, 1) for i, j in it.combinations(range(len(nodes)), 2)]
    # Nd, Cn, a, q = lo2(nodes, lines)
    # lo_plot(Nd, Cn, a, q, max(a) * 1e-3, str='Optimal layout', update=False, plane='xy')
    # print(a)  


