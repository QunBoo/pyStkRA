import numpy as np
import paraClass as para
import copy
import matplotlib.pyplot as plt

class optimizeHJTORA(para.optimizeAlgorithm):
    def __init__(self,para0):
        super(optimizeHJTORA, self).__init__(para0)

    def genOriginX(self):
        H = self.H
        [userNumber, serverNumber, sub_bandNumber] = H.shape
        seed = np.zeros([userNumber,serverNumber,sub_bandNumber])
        old_J = np.zeros([userNumber,serverNumber,sub_bandNumber])
        Pu_max_M = self.Pu_max * np.ones(userNumber)
        for user in range(0,userNumber):
            for server in range(0, serverNumber):
                for band in range(0, sub_bandNumber):
                    seed[user,server,band] = 1
                    du = self.Tu_data[user]
                    Ttol_buff = self.Ttol[user,server,band]
                    # T_com_Pu_max = self.Tcommu(seed, user,server,band)
                    Ttol_flag = self.check_Ttol(seed)
                    if(Ttol_flag == 1):
                        [old_J_temp, old_F, old_Pu] = self.RA(seed)
                        [old_J[user,server,band],old_F,old_Pu] = self.RA(seed)
                    else:
                        old_J[user,server,band] = 0
                    seed[user,server,band] = 0
        [user,server,band] = np.unravel_index(old_J.argmax(),old_J.shape)
        if(old_J[user,server,band] > 0):
            seed[user,server,band] = 1
        [J,F,Pu] = self.RA(seed)
        return seed,J,F,Pu

    def remove(self, X):
        [userNumber, serverNumber, sub_bandNumber] = X.shape
        user = 0
        server = 0
        band = 0
        not_find = 1
        [old_J, old_F, old_Pu] = self.RA(X)
        while(not_find == 1 & user != userNumber & server != serverNumber & band != sub_bandNumber):
            not_find = 1
            for user in range(0, userNumber):
                for server in range(0,serverNumber):
                    for band in range(0,sub_bandNumber):
                        if X[user,server,band] == 1:
                            X[user,server,band] = 0
                            [new_J, new_F, new_Pu] = self.RA(X)
                            if new_J > (1 + 0.001) * old_J:
                                not_find = 0
                                old_J = copy.deepcopy(new_J)
                                old_F = copy.deepcopy(new_F)
                                old_Pu = copy.deepcopy(new_Pu)
                                break
                            else:
                                X[user,server,band] = 1
                    if not_find == 0:
                        break
        return X, old_J, old_F, old_Pu, not_find

    def exchange(self,X):
        [userNumber, serverNumber, sub_bandNumber] = X.shape
        not_find = 1
        [old_J, old_F, old_Pu] = self.RA(X)
        X_new = copy.deepcopy(X)
        # X_check = X
        # X_check_flag1 = self.check_Ttol(X_check)
        for user in range(0,userNumber):
            for server in range(0,serverNumber):
                for band in range(0, sub_bandNumber):
                    if X[user,server,band] == 0:
                        X_new[user,...] = 0
                        X_new[...,server,band] = 0
                        X_new[user,server,band] = 1
                        [new_J,new_F,new_Pu] = self.RA(X_new)
                        # T_com_Pu_max = self.Tcommu(X_new,user,server,band,self.Pu_max)
                        Ttol_flag = self.check_Ttol(X_new)
                        # Ttol_buff = self.Ttol[user,server,band]
                        if (new_J > (1 + 0.001) * old_J) & (Ttol_flag == 1):
                            not_find = 0
                            old_J = copy.deepcopy(new_J)
                            old_F = copy.deepcopy(new_F)
                            old_Pu = copy.deepcopy(new_Pu)
                            X = copy.deepcopy(X_new)
                        else:
                            X_new = copy.deepcopy(X)
                if not_find == 0:
                    break
            if not_find == 0:
                break
        # X_check_flag2 = self.check_Ttol(X_check)
        # X_check_flag = (X_check == X).all()
        return X, old_J, old_F, old_Pu, not_find

    def ta(self):
        [X, J, F, Pu_out] = self.genOriginX()
        iterations = 1
        picture = np.array([0])
        flag = 1
        while(flag == 1):
            flag = 0
            [X, J, F, Pu_out, not_find_remove] = self.remove(X)
            if not_find_remove == 1:
                check1 = self.check_Ttol(X)
                [X, J, F, Pu_out, not_find_exchange] = self.exchange(X)
                check2 = self.check_Ttol(X)
                if(not_find_exchange == 0):
                    flag = 1
            # check = self.check_Ttol(X)
            iterations = iterations +1
            picture = np.append(picture,J)

        # pic_x = range(1,iterations)
        # pic_y = picture[1:]
        # plt.plot(pic_x, pic_y, linestyle='-', color='b',marker='^')
        # plt.show()
        # [J, F, Pu_out] = self.RA(X)
        return X, J, F, Pu_out



