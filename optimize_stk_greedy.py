import numpy as np
import paraClass as para

class optimizeGreedy(para.optimizeAlgorithm):
    def __init__(self,para0):
        super(optimizeGreedy, self).__init__(para0)

    def ta(self):
        H = self.H
        [userNumber, serverNumber, sub_bandNumber] = H.shape
        X = np.zeros([userNumber,serverNumber,sub_bandNumber])
        for user in range(0,userNumber):
            [_, server, _] = np.unravel_index(self.H.argmax(), self.H.shape)
            # np.argwhere(a == 0)
            sub_band = np.argwhere(X[:,server,:] == 0)
            if sub_band.size > 0:
                X[user,server,sub_band[0]] = 1
                check_Ttol_logi = self.check_Ttol(X)
                if check_Ttol_logi == 0:
                    X[user,server,sub_band[0]] = 0
        [J, F, Pu_out] = self.RA(X)
        return X, J, F, Pu_out

