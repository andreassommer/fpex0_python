import numpy as np
from fpex0.fokker_planck import FokkerPlanck
from scipy import sparse


class FokkerPlanckVDE(FokkerPlanck):

    @staticmethod
    def FokkerPlanckODE_VDE_inconsistent(self, t, Gp, u_sol, h, driftFcn, driftFcn_p, driftParams, diffusionFcn, diffusionFcn_p, diffusionParams, n_p, FPobj, IDobj):
        # t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams
        # in Setup.Parameters:
        # self.p0 = np.concatenate((p0_FPdrift, p0_FPdiffusion, p0_iniDist))

        # prepare variables
        np_FP = len(driftParams) + len(diffusionParams)
        np_IC = n_p - np_FP
        u  = u_sol(t) # call interpolated solution
        nu = len(u)
        Gp_matrix = Gp.reshape((nu,n_p), order='F') # convert Gp to matrix

        # determine VDE right-hand-side
        dfdp_FP = FPobj.FokkerPlanckODE_dfdp_FP(t, u, h, driftFcn_p, driftParams, diffusionFcn_p, diffusionParams) # RHS derivatives w.r.t. parameters
        dfdp = np.column_stack( (dfdp_FP, np.zeros((nu, np_IC))) )
        dfdu = self.FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams, banded=False)
        dGp  = (dfdu @ Gp_matrix).todense() + dfdp

        return dGp.flatten(order='F')
    

    def FokkerPlanckVDE_ODE(self, t, uu, h, driftFcn, driftFcn_p, driftParams, diffusionFcn, diffusionFcn_p, diffusionParams, n_p):        
        """
        Docs.
        """
        # prepare variables
        nuu = len(uu)
        np_FP = len(driftParams) + len(diffusionParams)
        np_IC = n_p - np_FP
        nu = nuu//(n_p+1)
        u  = uu[:nu]
        Gp = uu[nu:]
        Gp_matrix = Gp.reshape((nu,n_p), order='F') # convert Gp to matrix

        # nominal right-hand-side
        f = FokkerPlanck.FokkerPlanckODE(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams)
        
        # VDE right-hand-side
        dfdp_FP = self.FokkerPlanckODE_dfdp_FP(t, u, h, driftFcn_p, driftParams, diffusionFcn_p, diffusionParams) # RHS derivatives w.r.t. parameters
        dfdp = np.column_stack( (dfdp_FP, np.zeros((nu, np_IC))) )
        dfdu = self.FokkerPlanckODE_dfdu(t, u, h, driftFcn, driftParams, diffusionFcn, diffusionParams, banded=False)
        Gp_sparse = sparse.csr_matrix(Gp_matrix)
        dGp  = (dfdu @ Gp_sparse).todense() + dfdp

        # combine nominal and VDE right-hand-sides
        RHS = np.append(f, dGp.flatten(order='F'))
        return RHS


    @staticmethod
    def FokkerPlanckODE_VDE_dffduu():
        # not implemented. only needed if VDE should be solved implicitly.
        pass