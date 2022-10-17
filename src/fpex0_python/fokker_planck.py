import numpy as np
from scipy import sparse

class FokkerPlanck:
    def __init__(self):
        self._NN         = None # variables for FokkerPlanckODE_dfdu()
        self._A          = None
        self._B          = None

    @staticmethod
    def FokkerPlanckODE(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams):
        """
        ODE RHS of Fokker-Planck PDE by using the method of lines (MOL) approach.

        FP-PDE
        >       u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)

        FD-Approx
        >       u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
        >       u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h

        ## Takes           
        **t**: Time

        **x**: State vector
        
        **h**: MOL interval size
        
        **driftFcn**
        <br> Drift function, evaluated at t, driftP.
        
        **driftParams**
        <br> Parameter vector for drift function.
        
        **diffusionFcn**
        <br> Diffusion function, evaluated at t, driftP.
        
        **diffusionParams**
        <br> Parameter vector for diffusion function.

        ## Returns
        **dx**
        <br> rhs vector.
        """
    
        # number of grid points and number of simultaneously requested vectors
        shape = u.shape
        # shape is (N,) if u is a one-dim array
        if len(shape) > 1:
            (N, vectors) = shape
        else:
            N = shape[0]
            vectors = 1

        # preallocate vectors for A*u (approx. of u_x) and B*u (approx. of u_xx)
        Bu = np.zeros_like(u)
        Au = np.zeros_like(u)

        # numpy indexing is works differently for 1d and Nd arrays
        if vectors > 1:
            # first node
            Bu[0,:] = ( -2*u[0,:] + 2*u[1,:] ) / h**2
            # Au(1) is zero

            # inner nodes (remember: python slicing is exclusive right)
            i = np.arange(1,N-1)
            Au[i,:] = ( u[i-1,:] - u[i+1,:] ) / (2*h)            # 1st derivative stencil and scale
            Bu[i,:] = ( u[i-1,:] - 2*u[i,:] + u[i+1,:] ) / h**2  # 2nd derivative stencil and scale

            # last node
            Bu[N-1,:] = ( -2*u[N-1,:] + 2*u[N-2,:] ) / h**2
            # Au(N) is zero
        elif vectors == 1:
            # first node
            Bu[0] = ( -2*u[0] + 2*u[1] ) / h**2
            # Au(1) is zero

            # inner nodes (remember: python slicing is exclusive right)
            i = np.arange(1,N-1)
            Au[i] = ( u[i-1] - u[i+1] ) / (2*h)           # 1st derivative stencil and scale
            Bu[i] = ( u[i-1] - 2*u[i] + u[i+1] ) / h**2   # 2nd derivative stencil and scale

            # last node
            Bu[N-1] = ( -2*u[N-1] + 2*u[N-2] ) / h**2
            # Au[N-1] is zero      

        # evaluate drift and diffusion
        alpha = driftFcn(t,driftParams)
        D = diffusionFcn(t,diffusionParams)

        # assemble rhs
        dx = -alpha*Au + D*Bu
        return dx


    def FokkerPlanckODE_dfdu(self, t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams):
        """
        Jacobian of FokkerPlanckODE w.r.t. state u, i.e. df/du. 
        
        FP-PDE
        >       u_t = - a(t,p) * u_x(t,x)  +  D(t,p) * u_xx(t,x)
        >
        >       where a is the driftFcn, and b is the diffusionFcn
         
        FD-Approx  
        >       u_x  = ( u(t,x+h) - u(t,x-h ) / (2h)
        >       u_xx = ( u(t,x+h) - 2u(t,x) + u(t,x-h) ) / h
        
        Jacobian
        >       df/du = -a(t,p) * A + D(t,p) * B
        >
        >       where A is the first-derivative stencil [-1 0 1]
        >       and B is the second-derivative stencil [1 -2 1] plus Robin boundary
         
        ## Takes           
        **t**: Time
        
        **x**: State vector
        
        **h**: MOL interval size
        
        **driftFcn**
        <br> Drift function, evaluated at t, driftP.
        
        **driftParams**
        <br> Parameter vector for drift function.
        
        **diffusionFcn**
        <br> Diffusion function, evaluated at t, driftP.
        
        **diffusionParams**
        <br> Parameter vector for diffusion function.
        

        ## Returns
        **dfdu**
        <br> Sparse jacobian of FokkerPlanckODE (scipy.sparse).
        """
        # dimension
        N = len(u)

        # if dimension has changed, we have to reassemble the stencil
        if N != self._NN:
            e1 = np.ones(N)
            e0 = np.zeros(N)
            # 1st order stencil
            self._A =  sparse.dia_matrix(( np.array([e1, e0, -e1]), [-1, 0, 1] ), shape=(N, N))
            self._A = sparse.csr_matrix(self._A)   # change format to edit single entries
            self._A[0  ,1  ] = 0
            self._A[N-1,N-2] = 0
            # 2nd order stencil + robin boundary
            self._B = sparse.dia_matrix(( np.array([e1 , -2*e1 ,  e1]), [-1,0,1] ), shape=(N, N))
            self._B = sparse.csr_matrix(self._B)
            self._B[0  ,1  ] = 2
            self._B[N-1,N-2] = 2
            self._NN = N
        
        # evaluate drift and diffusion
        alpha = driftFcn(t, driftParams)
        D     = diffusionFcn(t, diffusionParams)

        # asseble the jacobian
        dfdu = -alpha * self._A / (2*h)  +  D * self._B / h**2
        return dfdu

    @staticmethod
    def defaultDriftFcn(t,p):
        """
        Default drift function used in FPEX0.

        ## Takes
        **t**
        <br> Current time / heating rate.

        **p**
        <br> Drift parameter vector (drift-parameters only!).


        ## Returns
        **fval**
        <br> Function value for drift.
        """

        # linear drift parametrization
        fval = p[0] + p[1] * t
        return fval
 
    @staticmethod
    def defaultDiffusionFcn(t,p,betamax):
        """
        Default diffusion function used in FPEX0.

        ## Takes
        **t**
        <br> Current time / heating rate.

        **p**
        <br> Diffusion parameter vector (diffusion-parameters only!).

        **betamax**
        <br> Maximum time / heat rate (used for ensuring non-negativity).


        ## Returns
        **fval**
        <br> Function value for diffusion.
        """

        # linear parametrization ensuring non-negativity for non-negative p1 and p2
        fval = p[0] + t * (p[1] - p[0]) / betamax

        return fval