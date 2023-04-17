import itertools
import numpy as np
from scipy import sparse
import time

class FokkerPlanck:
    def __init__(self):
        self._NN         = None # variables for FokkerPlanckODE_dfdu()
        self._A          = None
        self._B          = None

    def FD_stencils(self, u, h, recompute=False):
        N = len(u)

        # do the stencils have to be reassembled?
        # DEV/NOTE: For moving grids, the jacobian has to be recomputed regularly.
        if N != self._NN or recompute:
            
            # uniform grid
            if type(h) is float or type(h) is np.float64:
                # uniform grid
                e1 = np.ones(N)
                e0 = np.zeros(N)
                # 1st order stencil
                self._A = sparse.dia_matrix(( np.array([-e1/(2*h), e0/(2*h), e1/(2*h)]), [-1,0,1] ), shape=(N, N))
                self._A = sparse.csr_matrix(self._A)   # change format to edit single entries
                self._A[0,  1  ] = 0
                self._A[N-1,N-2] = 0
                # 2nd order stencil + robin boundary
                self._B = sparse.dia_matrix(( np.array([e1/h**2, -2*e1/h**2,  e1/h**2]), [-1,0,1] ), shape=(N, N))
                self._B = sparse.csr_matrix(self._B)
                self._B[0,  1  ] = 2
                self._B[N-1,N-2] = 2
                self._NN = N
            
            # non-uniform grid
            elif type(h) is np.ndarray:
                # non-uniform grid
                # help variables
                e1 = np.ones(N)
                h_ratio = h[1:]/h[:-1]

                # 1st order stencil
                scale = h[1:]*(1+h_ratio)
                dsub = np.concatenate((-h_ratio**2/scale,[0, 0]))
                d0   = np.concatenate(([0],-(1-h_ratio**2)/scale, [0]))
                dsup = np.concatenate(([0, 0], e1[:-2]/scale))
                self._A = sparse.dia_matrix(( np.array([dsub, d0, dsup]), [-1,0,1] ), shape=(N, N))
                
                # 2nd order stencil + robin boundary
                scale = h[1:]*h[:-1]*(1+h_ratio)
                dsub = np.concatenate((2*h_ratio/scale,[2/(h[-1]**2), 0]))
                d0   = np.concatenate(([-2/(h[0]**2)],-2*(1+h_ratio)/scale, [-2/(h[-1]**2)]))
                dsup = np.concatenate(([0, 2/(h[0]**2)], 2*e1[:-2]/scale))
                self._B = sparse.dia_matrix(( np.array([dsub, d0,  dsup]), [-1,0,1] ), shape=(N, N))
                self._NN = N

        return self._A, self._B
    

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

        # numpy indexing works differently for 1d and Nd arrays
        # TODO: find a more general way for this. add flag "uniform" for example.
        if type(h) is float or type(h) is np.float64:
            # uniform (spatial) grid
            if vectors > 1:
                # first node
                Bu[0,:] = ( -2*u[0,:] + 2*u[1,:] ) / h**2
                # Au(1) is zero

                # inner nodes 
                # (remember: python slicing is exclusive right)
                i = np.arange(1,N-1)
                Au[i,:] = ( u[i+1,:] - u[i-1,:] ) / (2*h)            # 1st derivative stencil and scale
                Bu[i,:] = ( u[i-1,:] - 2*u[i,:] + u[i+1,:] ) / h**2  # 2nd derivative stencil and scale

                # last node
                Bu[N-1,:] = ( -2*u[N-1,:] + 2*u[N-2,:] ) / h**2
                # Au(N) is zero
            elif vectors == 1:
                # first node
                Bu[0] = ( -2*u[0] + 2*u[1] ) / h**2
                # (Au(1) is zero)

                # inner nodes (remember: python slicing is exclusive right)
                i = np.arange(1,N-1)
                #Au[i] = ( u[i+1] - u[i-1] ) / (2*h)           # 1st derivative central stencil
                #Au[i] = ( u[i+1] - u[i] ) / (h)               # 1st derivative forward stencil
                Au[i] = ( u[i] - u[i-1] ) / (h)               # 1st derivative backward stencil
                Bu[i] = ( u[i-1] - 2*u[i] + u[i+1] ) / h**2   # 2nd derivative central stencil

                # last node
                Bu[N-1] = ( -2*u[N-1] + 2*u[N-2] ) / h**2
                # (Au[N-1] is zero)      

        elif type(h) is np.ndarray or type(h) is list:
            # non-uniform grid
            # DEV/NOTE we need to care better about cancellation! We have problems with that sometimes in experiments
            if vectors > 1:
                pass
            
            elif vectors == 1:
                # first node
                Bu[0] = ( -2*u[0] + 2*u[1] ) / h[0]**2
                # (Au(1) is zero)

                # inner nodes 
                # (difference schemes: Sundqvist and Veronis, 1970 (A simple finite-difference grid with non-constant intervals))
                i = np.arange(1,N-1)
                Au[i] = (u[i+1] - (h[i]/h[i-1])**2 * u[i-1] - (1-(h[i]/h[i-1])**2)*u[i] ) / (h[i]*(1+h[i]/h[i-1]))
                Bu[i] = 2*(u[i+1] + h[i]/h[i-1]*u[i-1] - (1+h[i]/h[i-1])*u[i])/(h[i]*h[i-1]*(1+h[i]/h[i-1]))
                
                # last node
                Bu[N-1] = ( -2*u[N-1] + 2*u[N-2] ) / h[-1]**2
                # (Au[N-1] is zero)      

        # evaluate drift and diffusion
        alpha = driftFcn(t,driftParams)
        D = diffusionFcn(t,diffusionParams)

        # assemble rhs
        dx = -alpha*Au + D*Bu
        return dx


    def FokkerPlanckODE_dfdu(self, t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams, banded=False):
        # TODO: Needs to be extended to uniform grids for the use of implicit integrators!
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
        
        **banded**
        <br> Should the jacobian be provided in banded form?  
        Format is described in [scipy.linalg.solve_banded](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_banded.html).

        ## Returns
        **dfdu**
        <br> Sparse jacobian of FokkerPlanckODE (scipy.sparse).
        """
        # obtain FD stencils of the discretization
        A, B = self.FD_stencils(u, h)
        
        # evaluate drift and diffusion
        alpha = driftFcn(t, driftParams)
        D     = diffusionFcn(t, diffusionParams)

        # assemble the jacobian
        dfdu = -alpha*A  +  D*B

        if banded:
            # dfdu is required in banded form
            dfdu_dia = sparse.dia_matrix(dfdu)
            # LSODA fortran code seems to require one additional row
            dfdu_banded = np.zeros((4,self._NN))
            data = dfdu_dia.data
            offsets = dfdu_dia.offsets

            # for n=3 diagonals, the banded form is (d_1, d_0, d_-1) 
            # where d_i are the diagonal arrays
            map = [1,0,-1]
            for i,j in itertools.product(range(3), range(3)):
                if offsets[i] == map[j]:
                    dfdu_banded[j,:] = data[i]
            return dfdu_banded

        return dfdu


    def FokkerPlanckODE_dfdp_FP(self, t, u, h, driftFcn_p, driftParams, diffusionFcn_p, diffusionParams):
        # dimensions
        n = len(u)
        np_FP = len(driftParams) + len(diffusionParams)

        # prepare variables
        # derivatives
        alpha_p = driftFcn_p(t, driftParams)
        D_p     = diffusionFcn_p(t, diffusionParams)

        # preallocate
        dfdp_FP = np.zeros((n, np_FP))

        # precompute
        A, B = self.FD_stencils(u, h)
        Au = A@u
        Bu = B@u
        
        # compute entries
        np_dr = len(driftParams)
        for i in range(np_dr):
            dfdp_FP[:,i] = -alpha_p[i]*Au

        for i in range(np_dr, np_FP):
            dfdp_FP[:,i] = D_p[i]*Bu

        return dfdp_FP


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