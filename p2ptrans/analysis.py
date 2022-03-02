from .config import *
from .utils import normal, find_uvw, rotate
from .display import printMatAndDir

def readCrystParam(crystfile):
    """ Reads the 4 possible parameters of the cystfile"""
    
    # Default values
    ccell1 = np.eye(3)
    ccell2 = np.eye(3)
    planehkl = [1,0,0]
    diruvw = [0,1,0]
    
    try:
        with open(crystfile,"r") as f:
            content = f.readlines()
    except FileNotFoundError:
            content = []

    for l in content:
        if l[0].rstrip() == "#":
            continue
        line = l.split('=')
        if len(line) == 2:
            if line[0].rstrip()=="ccell1":
                ccell1 = eval(line[1].rstrip())
            elif line[0].rstrip()=="ccell2":
                ccell2 = eval(line[1].rstrip())
            elif line[0].rstrip()=="planehkl":
                planehkl = eval(line[1].rstrip())
            elif line[0].rstrip()=="diruvw":
                diruvw = eval(line[1].rstrip())
            else:
                print("WARNING: %s is not a supported input"%(line[0].rstrip()))
        elif len(line) > 2:
            raise SyntaxError(l)

    return ccell1, ccell2, planehkl, diruvw

def strainDirs(tmat, ftf=True):

    if ftf:
        eigval, P = la.eigh(tmat.T.dot(tmat))
        eigval = np.sqrt(eigval)
        
        idx = np.argsort(eigval)
        eigval = eigval[idx]
        P = P[:,idx]
        
        U = P.dot(np.diag(eigval)).dot(P.T)

        invEigval, Q = la.eigh(la.inv(tmat).T.dot(la.inv(tmat)))
        invEigval = np.sqrt(invEigval)
        idx = np.argsort(1/invEigval)
        Q = Q[:,idx]
        
    else:
        U = rotate(tmat, np.eye(3)).dot(tmat)

        eigval, P = la.eig(U)

        idx = np.argsort(eigval)
        eigval = eigval[idx]
        P = P[:,idx]
        
        iU = rotate(la.inv(tmat), np.eye(3)).dot(la.inv(tmat))
        
        invEigval, Q = la.eig(iU)
        invEigval = np.sqrt(invEigval)
        idx = np.argsort(1/invEigval)
        Q = Q[:,idx]
        
    P = normal(P)
    
    Q = normal(Q)

    return eigval, U, P, Q

def findHabit(U, P, eigval):
    
    # Using uniformly strained plane
    if abs(eigval[1]**2 - eigval[0]**2) > tol: 
        ratio = np.sqrt(abs((eigval[2]**2 - eigval[1]**2)/(eigval[1]**2 - eigval[0]**2)))
    
        # uvw direction of normal in cartesian coord
        planeHab = np.zeros((3,2))
        planeHab[:,0] = P[:,0] + ratio*P[:,2]
        planeHab[:,1] = P[:,0] - ratio*P[:,2]

    elif abs((eigval[2]**2 - eigval[1]**2)) < tol:
        # Any plane is good because all three strains are equal
        planeHab = np.zeros((3,2))
        planeHab[:,0] = P[:,0]
        planeHab[:,1] = -P[:,0]
        ratio = 1
    else:
        planeHab = np.zeros((3,2))
        planeHab[:,0] = P[:,2]
        planeHab[:,1] = -P[:,2]
        ratio = np.inf
    
    return planeHab, ratio

def findR(U, P=None, planeHab=None, ratio=None):

    if np.all(U==np.eye(3)):   
        return np.array([U,U])
    
    if planeHab is None or ratio is None or P is None:
        eigval, U, P, Q = strainDirs(U)
        planeHab, ratio = findHabit(U, P, eigval)
    
    V = np.zeros((2,3,3))
    M = np.zeros((2,3,3))
    R = np.zeros((2,3,3))
    
    for i in range(2):
        if np.isinf(ratio):
            V[i,:,0] = P[:,0]
        else:
            V[i,:,0] = ratio*P[:,0] - (1-2*i)*P[:,2]
            V[i,:,0] = V[i,:,0]/la.norm(V[i,:,0])
        V[i,:,1] = -(1-2*i)*P[:,1]
        V[i,:,2] = planeHab[:,i]/la.norm(planeHab[:,i])
        
        M[i,:,:2] = U.dot(V[i,:,:2])
        M[i,:,:2] = M[i,:,:2]/la.norm(M[i,:,0])
        M[i,:,2] = np.cross(M[i,:,0], M[i,:,1])
        M[i,:,2] = M[i,:,2] / la.norm(M[i,:,2])

        R[i,:,:] = V[i,:,:].dot(M[i,:,:].T)
        
    return R

def crystallography(tmat, A, B, ccellA=np.eye(3), ccellB=np.eye(3), planehkl=[1,0,0], diruvw=[1,0,0], fileA="input 1", fileB="input 2", ftf=True):

    """ This function does the crystallographic analysis using the deformation matrix tmat. It displays 
    strain directions, habit plane and OR given that tmat.dot(TC^(A)) = TC^(B). If you are using the result
    from findMatching and you want the strain directions from the initial to the final structure you must use
    la.inv(tmat) and inverse the order of A and B. The idea is that tmat, the input, is the transformation
    gradient matrix between A and B"""
    
    print("----------CRYSTALLOGRAPHY----------")
    print()
    eigval, U, P, Q = strainDirs(tmat, ftf=ftf)

    print("Strain Directions in %s (%s) coordinates:"%(B.name, fileB))
    printMatAndDir(P, ccellB)
    print()
    print("Strain Directions in %s (%s) coordinates:"%(A.name, fileA))
    printMatAndDir(Q, ccellA)
    print()
    print("Strains + 1 (eigenvalues)")
    print("    e1    e2    e3    ")
    print(' '.join(["% 5.3f"%(val) for val in eigval])) 
    print()

    planeHab, ratio = findHabit(U, P, eigval)
        
    print("Uniformly strained planes:")
    print()
    print("Exact plane hkl in %s (%s) coordinates:"%(B.name, fileB))
    print("(+): (% 6.4f, % 6.4f, % 6.4f)"%(*ccellB.T.dot(planeHab[:,0]),))
    print("(-): (% 6.4f, % 6.4f, % 6.4f)"%(*ccellB.T.dot(planeHab[:,1]),))
    print()
    print("Closest hkl:")
    print("(+): (% d, % d, % d)"%(*find_uvw(planeHab[:,0:1], la.inv(ccellB.T)),))
    print("(-): (% d, % d, % d)"%(*find_uvw(planeHab[:,1:2], la.inv(ccellB.T)),))
    print()

    R = findR(U, P=P, planeHab=planeHab, ratio=ratio)
    
    print("Orientation Relationship with habit plane:")
    print()
    orMat = np.zeros((2,3,3))
    resPlanehkl = np.zeros((2,3))
    resDiruvw = np.zeros((2,3))
    for i in range(2):
        orMat[i,:,:] = Q.dot(P.T).dot(R[i,:,:].T)
        resPlanehkl[i,:] = ccellB.T.dot(orMat[i,:,:].dot(la.inv(ccellA.T).dot(planehkl)))
        resDiruvw[i,:] = la.inv(ccellB).dot(orMat[i,:,:].dot(ccellA.dot(diruvw)))
    print("%s (%s) // %s (%s)"%(B.name, fileB, A.name, fileA))
    print("(+): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl[0,:], *resDiruvw[0,:]))
    print("(-): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl[1,:], *resDiruvw[1,:]))
    print()
    print("Approximate low index OR")
    resPlanehklClose = np.zeros((2,3))
    resDiruvwClose = np.zeros((2,3))
    for i in range(2):
        resPlanehklClose[i,:] = find_uvw(la.inv(ccellB.T).dot(resPlanehkl[i,:].reshape((3,1))), la.inv(ccellB.T)).T[0]
        resDiruvwClose[i,:] = find_uvw(ccellB.dot(resDiruvw[i,:].reshape((3,1))), ccellB).T[0]
    print("(+): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose[0,:], *resDiruvwClose[0,:]))
    print("(-): (% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose[1,:], *resDiruvwClose[1,:]))
    print()
    
    print("Orientation Relationship in thin films:")
    print()
    orMat = Q.dot(P.T)
    resPlanehkl = ccellB.T.dot(orMat.dot(la.inv(ccellA.T).dot(planehkl)))
    resDiruvw = la.inv(ccellB).dot(orMat.dot(ccellA.dot(diruvw)))
    print("%s (%s) // %s (%s)"%(B.name, fileB, A.name, fileA))
    print("(% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] // (% 6.4f, % 6.4f, % 6.4f) [% 6.4f, % 6.4f, % 6.4f]"%(*planehkl, *diruvw, *resPlanehkl, *resDiruvw))
    print()
    print("Approximate low index OR")
    resPlanehklClose = find_uvw(la.inv(ccellB.T).dot(resPlanehkl.reshape((3,1))), la.inv(ccellB.T)).T[0]
    resDiruvwClose = find_uvw(ccellB.dot(resDiruvw.reshape((3,1))), ccellB).T[0]
    print("(% 2d, % 2d, % 2d) [% 2d, % 2d, % 2d] //  (% d, % d, % d) [% d, % d, % d]"%(*planehkl, *diruvw, *resPlanehklClose, *resDiruvwClose))
    print()

    return eigval, U, P, Q, planeHab 

