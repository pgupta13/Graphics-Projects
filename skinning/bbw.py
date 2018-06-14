from numpy import *
import scipy.sparse.linalg
#import numpy as np

def cvxopt_sparse( scipy_sparse_matrix ):
    '''
    Converts a scipy.sparse matrix to a cvxopt sparse matrix.
    '''
    import cvxopt
    P = scipy.sparse.coo_matrix( scipy_sparse_matrix )
    M = cvxopt.spmatrix( P.data, asarray( P.row, dtype = int ), asarray( P.col, dtype = int ), P.shape )
    return M

def linear_blend_skin_2D( vertices, weights, transforms ):

    '''
    Given an array of 2D `vertices`,
    an array of #vertices-by-#transforms `weights`,
    and an array of 3x3 matrices `transforms`,
    returns the linear blend skinning deformation of each vertex.
    
    NOTE: The input 2D vertices are not homogeneous.
          They are just [ [ x0, y0 ], [ x1, y1 ], ... ].
          You must convert them to homogeneous coordinates to see the effects of translation.
          You must return non-homogeneous 2D coordinates.
    '''
    
    ## Add your code here.

    x = len(vertices)
    #print(str(x)+'vert')
    #on= ones((x,1))
    #vert = append(vertices,on,axis=1)
    vert = vertices
    vert[:,2]=1;
    #print(vert)

    A = [0] * x
    for i in range(x):
        A[i] = [0] * 2
    #print(A)

    y = len(weights[0])
    for i in range(x):
        temp = zeros((3,3))
        for j in range(y):
            temp2 = weights[i,j]*transforms[j]
            temp=temp+temp2
        A[i]=temp@(vertices[i][0],vertices[i][1],1.)

    A=delete(A,2,1)
    #print(A.shape)
    return array( A )

def laplacian_and_mass_matrices( faces, vertices ):
    '''
    Given `vertices`, an array of N 2D points pi = [xi, yi] (equivalently, an N-by-2 array),
    and `faces`, an array of F triplets of integer indices into vertices,
    where the triplet faces[f][0], faces[f][1], faces[f][2]
    are the indices of the three vertices that make up triangle f,
    return two N-by-N sparse matrices: Laplacian, Mass.
    '''

    N = len(vertices)



    #A = [0] * N
    #for i in range(N):
    #    A[i] = [0] * N

    A= zeros((N,N))
    #print(A)

    F=len(faces)
    for f in range(F):
        Ind1=faces[f][0]
        Ind2=faces[f][1]
        Ind3=faces[f][2]

        A[Ind1][Ind2]=1; A[Ind2][Ind1]=1;
        A[Ind1][Ind3]=1;A[Ind3][Ind1]=1;
        A[Ind2][Ind3]=1;A[Ind3][Ind2]=1;

    #M = [0] * N
    #for i in range(N):
    #    M[i] = [0] * N
    M= zeros((N,N))

    for i in range(N):
        M[i][i]=sum(A[i])

    #print(M)

    subtractRes = array(M) - array(A)
    sparOfM= scipy.sparse.csc_matrix(M)
    invOfM = scipy.sparse.linalg.inv(sparOfM)
    Lap = scipy.sparse.csc_matrix(invOfM@subtractRes)
    #sparceM= scipy.sparse.csc_matrix(M)
    #print(Lap)

    return Lap, sparOfM

def bbw( faces, vertices, handles, laplacian_mode, solver_mode ):
    '''
    Given `faces` and `vertices` as would be passed to laplacian_and_mass_matrices(),
    an array of H integer indices into `vertices` representing the handle vertices,
    a string `laplacian_mode` which will be one of "graph" or "cotangent",
    and a string `solver_mode` which will be one of "bounded" or "unbounded",
    return an #vertices-by-#handles weight matrix W, where W[i,j] is the influence weight
    of the j-th handle on the i-th vertex.
    Each row of W must sum to 1.
    If mode == "bounded", apply inequality constraints so that the weights are
    all between 0 and 1.
    '''
    

    import cvxopt
    L,M = laplacian_and_mass_matrices(faces,vertices)
    P = L.T@M@L             #B N x N

    #P = cvxopt.matrix(P,tc='d')
    P= cvxopt_sparse(P)



    N = len( vertices )
    H = len( handles )

    q = zeros((N,1)) #q - N x 1
    q = cvxopt.matrix(q,tc='d')
    Aieq=zeros((2*N,N))
    Aieq[0:N,0:N]= identity(N) #G - 2N x N
    Aieq[N:2*N,0:N]= -1*identity(N)
    Aieq = cvxopt.matrix(Aieq,tc='d')
    #Aieq=cvxopt_sparse(Aieq)
    Bieq=zeros((2*N,1))
    Bieq[0:N,0]=ones(N) #h 2N x 1
    Bieq=cvxopt.matrix(Bieq,tc='d')
    #Bieq=cvxopt_sparse(Bieq)
    beq=zeros((H,1)) #b H x 1
    #beq=cvxopt_sparse(beq)

    Aeq=zeros((H,N)) #A H x N



    for i in range(H):
        Aeq[i][handles[i]] = 1
    #Aeq = cvxopt.matrix(Aeq, tc='d')
    Aeq=cvxopt_sparse(Aeq)

    weights= ones((N,H))

    if(solver_mode=="bounded"):
        for i in range(H):
            beq[i] = 1
            # P, q, G, h, A, b
            z=cvxopt.matrix(beq,tc='d')
            result = cvxopt.solvers.qp(P, q, G=Aieq, h=Bieq, A = Aeq, b=z)
            weights[:, i] = array(result['x'])[:, 0]
            beq[i]=0
    else:
        for i in range(H):
            beq[i] = 1
            # P, q, G, h, A, b
            z = cvxopt.matrix(beq, tc='d')
            result = cvxopt.solvers.qp(P, q, A=Aeq, b=z)
            weights[:, i] = array(result['x'])[:, 0]
            beq[i] = 0


    rowSum = weights.sum(axis=1)
    weights_new=weights/rowSum[:,newaxis]
    #weights = (1./H)*weights
    print(weights_new)

    return weights_new
