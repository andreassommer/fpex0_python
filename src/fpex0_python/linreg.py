import numpy as np


class LinearRegression():
    """
    Assumes a linear model:  Yi = a + b*Xi + ei    with ei ~ N(0,s^2)
    and calculates the simple linear regression according to [1].

    Executes regression and holds regression parameters. 
    <br> If a first regression has been made, any following execution of self.linReg(X, Y) will update the
    regression with the additional point (X, Y).
    
    ## Parameters
    Regression represented by following fields:
    **a** Estimate for a
    
    **b**: Estimate for b
    
    **s2**: Estimate for s^2
    
    **Xmean**: Mean value of Xi
    
    **Ymean**: mean value of Yi
    
    **Sxx**: S_X^2
    
    **Sxy**: S_XY
    
    **Syy**: S_Y^2
    
    **See**: S_e^2


    ### Sources
    [1] Jerome H. Klotz: "Updating Simple Linear Regression", 1995. <br>
    [2] Statistica Sinica 5 (1995), 399-403.
    
    
    """
    def __init__(self):
        self.n     = 0
        self.Xmean = None
        self.Ymean = None
        self.See   = None
        self.Sxx   = None
        self.Sxy   = None
        self.a     = None
        self.b     = None
        self.s2    = None


    def linReg(self, X, Y):
        """
        Sets or updates the regression attributes.

        ## Takes
        **X**: numpy.array or number <br>
        **Y**: numpy.array or number

        X, Y must be either both arrays or both numbers.

        ### Note
        A vector update has not yet been implemented, but is planned.
        """


        if self.n == 0:
            n = len(X)
            Xmean = np.mean(X)
            Ymean = np.mean(Y)
            XminusXmean = X - Xmean
            YminusYmean = Y - Ymean
            Sxx = np.dot(XminusXmean,XminusXmean)
            Sxy = np.dot(XminusXmean,YminusYmean)
            b   = Sxy / Sxx
            a   = Ymean - b*Xmean
            See = YminusYmean - b*XminusXmean
            See = np.dot(See, See)
            s2  = See / (n-2)
            
            # assign values
            self.n     = n
            self.Xmean = Xmean
            self.Ymean = Ymean
            self.See   = See
            self.Sxx   = Sxx
            self.Sxy   = Sxy
            self.a     = a
            self.b     = b
            self.s2    = s2


        # is X primitive?  
        # --> this also accepts strings for example,
        #     but gives the user the freedom to use customized floats
        elif not hasattr(X, '__len__'):
            # accessors
            n      = self.n
            bn     = self.b
            Xmeann = self.Xmean
            Ymeann = self.Ymean
            Sxxn   = self.Sxx
            Sxyn   = self.Sxy
            Seen   = self.See
    

            ## scalar update
            N = n + 1

            XmeanN = Xmeann + (X - Xmeann)/N
            YmeanN = Ymeann + (Y - Ymeann)/N

            XminusXmean = X - Xmeann
            YminusYmean = Y - Ymeann

            SxxN = Sxxn + n/N * XminusXmean * XminusXmean
            SxyN = Sxyn + n/N * XminusXmean * YminusYmean
            SeeN = Seen + n/N * (YminusYmean - bn*XminusXmean)**2 * Sxxn / SxxN

            bN = SxyN / SxxN
            aN = YmeanN - bN*XmeanN  # Formula (1)
            s2N = SeeN / (N-2)

            ## store values
            self.n     = N
            self.Xmean = XmeanN
            self.Ymean = YmeanN
            self.Sxx   = SxxN
            self.Sxy   = SxyN
            self.See   = SeeN
            self.a     = aN
            self.b     = bN
            self.s2    = s2N


        elif len(X) >  1:
            # accessors
            n      = self.n
            bn     = self.b
            Xmeann = self.Xmean
            Ymeann = self.Ymean
            Sxxn   = self.Sxx
            Sxyn   = self.Sxy
            Seen   = self.See
    
            nnew = len(X)
            N = n + nnew
            sqrtn = np.sqrt(n)
            sqrtN = np.sqrt(N)

            XmeanN = (n / N) * (Xmeann + nnew/n * np.mean(X))  # own derivation
            YmeanN = (n / N) * (Ymeann + nnew/n * np.mean(Y))  #
         
            XmeanSTAR = (Xmeann*sqrtn + XmeanN*sqrtN) / (sqrtn + sqrtN)  # below Formula (6)
            YmeanSTAR = (Ymeann*sqrtn + YmeanN*sqrtN) / (sqrtn + sqrtN)  # below Formula (7)
         
            XminusXmeanSTAR = (X - XmeanSTAR)
            SxxSTAR = np.dot(XminusXmeanSTAR, XminusXmeanSTAR)
            SxxN    = Sxxn + SxxSTAR                            # Formula (6)
         
            YminusYmeanSTAR = (Y - YmeanSTAR)                   #
            SxySTAR = np.dot(XminusXmeanSTAR, YminusYmeanSTAR)         #
            SxyN    = Sxyn + SxySTAR                            # Formula (7)
         
            bN = SxyN / SxxN                                    # Formula (8)
         
            Sxn = np.sqrt(Sxxn)
            SxN = np.sqrt(SxxN)
            bSTAR = (bn * Sxn + bN*SxN) / (Sxn + SxN)           # Formula (10)
         
            SeeSTAR = YminusYmeanSTAR - bSTAR*XminusXmeanSTAR   #
            SeeSTAR = np.dot(SeeSTAR, SeeSTAR)                         # Formula (9)
            SeeN    = Seen + SeeSTAR                            # Formula (8)
         
            aN  = YmeanN - bN*XmeanN   # Formula (1)
            s2N = SeeN / (N-2)


            ## store values
            self.n     = N
            self.Xmean = XmeanN
            self.Ymean = YmeanN
            self.Sxx   = SxxN
            self.Sxy   = SxyN
            self.See   = SeeN
            self.a     = aN
            self.b     = bN
            self.s2    = s2N

