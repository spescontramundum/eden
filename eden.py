
import math
import numpy as np
import random as rand
from matplotlib import pyplot as plt
from matplotlib import patches
#
# # The following is an instance of the Eden growth model for a two dimensional lattice
# # The process iteratively constructs a set with a subset called the frontier
# # that grows under simple stochastic dynamics.  It each time, for each element
# # of the frontier, nearest neighbors are added independently with fixed probability
# # This instantiation defines the frontier to be any point in the set with a
# # nearest neighbor not in the set.

prob = .2

#generates points in each quadrant
def quad_set(i,j):
    return [(i, j), (-i, j),(-i, -j),(i, -j)]


#create a dictionary to store the values of points
def create_base(n):
    M = {}

    for i in range(n):
        for j in range(n):
            for k in quad_set(i, j):
                M[k] = 0
    return M


def NearPoint(i, j, M):
    r = np.NAN
    t = np.NAN
    l = np.NAN
    b = np.NAN

    # this seemed easier than doing a bunch of edge cases
    try:
        # value of point to right
        r = M[(i + 1, j)]
    except:
        pass
    try:
        # value of point above
        t = M[(i, j + 1)]
    except:
        pass
    try:
        # value of point to left
        l = M[(i - 1, j)]
    except:
        pass
    try:
        # value of point below
        b = M[(i, j - 1)]
    except:
        pass

    return r, t, l, b



def UpdatePoint(point_list):
    ''' update any point with neighboring points in the set, making it part of the set
    with a random probability'''
    r = point_list[0]
    t = point_list[1]
    l = point_list[2]
    b = point_list[3]

    def roll():
        return rand.random()

    if r == 1 and roll() < prob:
        return 1
    if t == 1 and roll() < prob:
        return 1
    if l == 1 and roll() < prob:
        return 1
    if b == 1 and roll() < prob:
        return 1

    # if none of the rolls are successful return 0
    return 0

def UpdateModel(M,n):

    for i in range(n):
        for j in range(n):
            for k in quad_set(i,j):
                p = M[k]
                a = list(k)[0]
                b = list(k)[1]
                if p == 0:
                    M[k] = UpdatePoint(list(NearPoint(a,b, M)))

    return M


def RunModel(n):
    '''creates a set M which updates (n times) the exterior points using UpdateModel() and UpdatePoint()'''
    M = create_base(n)
    #random starting seed
    M[(0,0)] = 1
    for i in range(n):
        M = UpdateModel(M,n)
    return M



def find_edges(n):
    M = RunModel(n)
    interior = {}
    edge_set = {}

    for i in range(n):
        for j in range(n):
            for k in [(i,j),(-i, j),(-i, -j), (i, -j)]:
                p = M[k]
                a = list(k)[0]
                b = list(k)[1]
                if p == 1:
                    point_list = list(NearPoint(a, b, M))
                    if np.nansum(point_list) == 4:
                        interior[k] = M[k]
                    else:
                        edge_set[k] = M[k]

    return edge_set, interior






A, B = find_edges(14)




def check_rin(C):
    ''' checks within a circle, returns False if there are any points not belonging to the Eden set'''
    for c in C:
        if C[c] == 0:
            return False
    return True

def circle_set(M, r):
    C = {}
    for i in range(int(r) +1):
        for j in range(int(r) +1):
            for k in quad_set(i,j):
                a = list(k)[0]
                b = list(k)[1]
                if a**2 + b**2 <= r**2:
                    try:
                        C[k] = M[k]
                    except:
                        pass
    return C


def Rin(M,n, precision):
    '''calculates inner radius'''
    r = n
    r_old = None
    for i in range(1, precision):
        C = circle_set(M, r)
        r_old = r
        if check_rin(C) == True:
            # if it includes only Eden set points try a bigger radius
            r = (r + n / (2** i))
        else:
            # if the circle contains points not in the Eden set try a smaller one
            r = (r - n / (2 ** i))

    if check_rin(circle_set(M, r)) == False:
        return r_old
    else:
        return r



M = RunModel(14)
r_in = Rin(M, 14, 8)
if __name__=="__main__":

    Ex = []
    Ey = []
    Ix = []
    Iy = []

    Akeys = list(A.keys())
    Bkeys = list(B.keys())

    for a in range(len(A)):
        la = list(Akeys[a])
        Ex.append(la[0])
        Ey.append(la[1])

    for b in range(len(B)):
        lb = list(Bkeys[b])
        Ix.append(lb[0])
        Iy.append(lb[1])

    mass = len(Ex)
    cgx = np.sum(Ex)/mass
    cgy = np.sum(Ey)/mass
    print(cgx,cgy)

    plt.scatter(Ex,Ey, s=.5, c='red')
    plt.scatter(Ix,Iy, s=.5, c='blue')
    c_in = plt.Circle((0, 0), radius=r_in, edgecolor="blue", fill=False)
    fig = plt.gcf()
    ax = fig.gca()
    ax.add_patch(c_in)
    plt.title("Eden Growth Model")
    plt.figure(figsize=(20, 20))
    plt.show()


