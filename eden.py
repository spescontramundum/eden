
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

prob = .1

#generates points in each quadrant
def quad_set(i,j):
    return [(i, j), (-i, j),(-i, -j),(i, -j)]


#create a dictionary to store the values of points
def create_base(n):
    M = {}
    bounds = 4*n #int(np.sqrt(2)* n) + 10
    for i in range(bounds):
        for j in range(bounds):
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
    bounds = 3*n #int(np.sqrt(2)* n + 1)
    M_new = M
    for i in range(bounds):
        for j in range(bounds):
            for k in quad_set(i,j):
                p = M[k]
                a = list(k)[0]
                b = list(k)[1]
                if p == 0:
                    M_new[k] = UpdatePoint(list(NearPoint(a,b, M)))

    return M_new


def RunModel(n):
    '''creates a set M which updates (n times) the exterior points using UpdateModel() and UpdatePoint()'''
    M = create_base(n)
    #random starting seed
    M[(0,0)] = 1
    for i in range(n):
        M = UpdateModel(M,n)
    return M



def find_edges(M,n):

    interior = {}
    edge_set = {}
    bounds = 3*n #int(np.sqrt(2) * n) + 1
    for i in range(bounds):
        for j in range(bounds):
            for k in quad_set(i,j):
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


def check_rin(C):
    ''' checks within a circle, returns False if there are any points not belonging to the Eden set'''
    for c in C:
        if C[c] == 0:
            return False
    return True

def circle_set(M, r):
    C = {}
    for i in range(int(np.sqrt(2)*r)+1):
        for j in range(int(np.sqrt(2)*r)+1):
            for k in quad_set(i,j):
                a = list(k)[0]
                b = list(k)[1]
                #if (a - cg[0])**2 + (b - cg[1])**2 <= r**2:
                if a**2 + b**2 <= r**2:
                    C[k] = M[k]


    return C


def Rin(M, n, precision=6):
    '''calculates inner radius'''
    R = int(np.sqrt(2)*n + 1)
    r = R
    r_list = []
    for i in range(1, precision):

        C = circle_set(M, r)
        r_list.append(r)

        if check_rin(C) == True:
            # if it includes only Eden set points try a bigger radius
            r = (r + R / (2 ** i))
        else:
            # if the circle contains points not in the Eden set try a smaller one
            r = (r - R / (2 ** i))

    for j in range(1, precision):
        if check_rin(circle_set(M, r_list[-j])) == True:
            return r_list[-j]


def Rout(M, n, precision=5):
    '''calculates inner radius'''
    #set starting radius
    R = int(np.sqrt(2)*n + 1)
    m = np.sum(list(M.values()))

    r_list = [R]
    includes = False


    for i in range(1, precision):
        r = r_list[-1]
        C = circle_set(M, r)

        count = 0
        for c in C:
            if C[c] == 1:
                count += 1

        if count == m:
            # if the circle does not contain all the points try a larger one
            r = (r - (R / (2 ** i)))
            r_list.append(r)
            includes = True

        else:
            r = (r + (R / (2 ** i)))
            r_list.append(r)
            includes = False

    if includes == False:
        return r_list[-1]
    else:
        return r


def estimate_r(n, trials = 10):

    r_ins = []
    r_outs = []
    for i in range(trials):
        M = RunModel(n)
        r_ins.append(Rin(M, n, precision=12))
        r_outs.append(Rout(M, n, precision=12))
    rin_avg = np.sum(r_ins)/trials
    rout_avg = np.sum(r_outs)/trials
    return rin_avg,rout_avg


def plot_radii(n,r_ins, r_outs):

    plt.plot(range(n), r_ins, color='light')
    plt.plot(range(n), r_outs, color='red')
    plt.title("In-radii(b) and out-radii(r) as a function of growth epochs")
    plt.figure(figsize=(20, 20))
    plt.show()

    #print("The estimate for the  in-radius of a model with ", n, " steps is: ", rin_avg )

def plot_eccent(n, E_t):
    plt.plot(range(n), E_t, color='purple')
    plt.title("In-radii(b) and out-radii(r) as a function of growth epochs")
    plt.figure(figsize=(20, 20))
    plt.show()

def long_term(n):
    r_ins = []
    r_outs = []

    for i in range(n):
        est = list(estimate_r(n))
        r_ins.append(est[0])
        r_outs.append(est[1])

    plot_radii(n, r_ins,r_outs)

    E_t = []
    for i in range(n):
        E_t.append(r_ins[i]/r_outs[i])
    plot_eccent(n,E_t)


def plot_Model(n):

    M = RunModel(n)
    A, B = find_edges(M, n)


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

    # I don't know whether the center of mass helps, and haven't run the numbers yet
    # my guess is probably not, but it was the only potential easy fix I could think of
    # for in radii that left a bunch of stuff out
    # mass = len(Ex)
    # cgx = np.sum(Ex)/mass
    # cgy = np.sum(Ey)/mass

    r_in = Rin(M, n, precision=12)
    r_out = Rout(M, n, precision=12)




    plt.scatter(Ex,Ey, s=.5, c='red')
    plt.scatter(Ix,Iy, s=.5, c='blue')
    c_in = plt.Circle((0,0), radius=r_in, edgecolor="midnightblue", fill=False)
    c_out = plt.Circle((0, 0), radius=r_out, edgecolor="red", fill=False)
    fig = plt.gcf()
    ax = fig.gca()
    fig2 = plt.gcf()
    ax2 = fig.gca()
    ax.add_patch(c_in)
    ax2.add_patch(c_out)
    plt.title("Eden Growth Model")
    plt.figure(figsize=(20, 20))
    plt.show()

if __name__=="__main__":
    plot_Model(14)
    #long_term(25)


