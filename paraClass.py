import numpy as np
import math
import scipy.io as io


class paraClass0:

    def __init__(self, beta_time, beta_enengy, Tu_data, Tu_circle,
                 tu_local, Eu_local, W, H, lamda, Pu, Sigma_square, Fs,
                 Eta_user, Pu_max, Pu_min, Ttol):
        self.beta_time = beta_time
        self.beta_enengy = beta_enengy
        self.Tu_data = Tu_data
        self.Tu_circle = Tu_circle
        self.tu_local = tu_local
        self.Eu_local = Eu_local
        self.W = W
        self.H = H
        self.lamda = lamda
        self.Pu = Pu
        self.Sigma_square = Sigma_square
        self.Fs = Fs
        self.Eta_user = Eta_user
        self.Pu_max = Pu_max
        self.Pu_min = Pu_min
        self.Ttol = Ttol


# classoptimizaAlgorithm
class optimizeAlgorithm(object):

    def __init__(self, para0):
        self.beta_time = para0.beta_time
        self.beta_enengy = para0.beta_enengy
        self.Tu_data = para0.Tu_data
        self.Tu_circle = para0.Tu_circle
        self.tu_local = para0.tu_local
        self.Eu_local = para0.Eu_local
        self.W = para0.W
        self.H = para0.H
        self.lamda = para0.lamda
        self.Pu = para0.Pu
        self.Sigma_square = para0.Sigma_square
        self.Fs = para0.Fs
        self.Eta_user = para0.Eta_user
        self.Pu_max = para0.Pu_max
        self.Pu_min = para0.Pu_min
        self.Ttol = para0.Ttol

    def genUs(self, x, server):
        #     生成对应服务器的用户矩阵
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        num = 0
        Us = np.array([])
        for user in range(0, userNumber):
            for sub_band in range(0, sub_bandNumber):
                if (x[user, server, sub_band] > 0):
                    num = num + 1
                    Us_temp = np.array([user, sub_band])
                    if (num == 1):
                        Us = Us_temp
                    else:
                        Us = np.row_stack((Us, Us_temp))

        return Us, num

    def cra(self, x):
        # CRA计算资源分配函数
        # x应为np.ndarray变量
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        Eta_user = self.Eta_user
        Fs = self.Fs
        F = np.zeros([userNumber, serverNumber])
        T = 0
        for server in range(0, serverNumber):
            [Us, n] = self.genUs(x, server)
            if n > 0:
                EtaRoot_sum = 0
                for i in range(0, n):
                    if n == 1:
                        eta_user_ptemp = Us[0]
                    else:
                        eta_user_ptemp = Us[i, 0]
                    eta_user_temp = Eta_user[eta_user_ptemp]
                    EtaRoot_sum = EtaRoot_sum + eta_user_temp ** 0.5
                for i in range(0, n):
                    if n == 1:
                        eta_user_ptemp = Us[0]
                    else:
                        eta_user_ptemp = Us[i, 0]
                    eta_user_temp = Eta_user[eta_user_ptemp]
                    F[eta_user_ptemp, server] = Fs[server] * eta_user_temp ** 0.5 / EtaRoot_sum
                T = T + 1 / Fs[server] * EtaRoot_sum ** 2
        res_cra = T
        return F, res_cra

    def getGamma(self, x, user, server, band, pu = float('nan')):
        # %GetGamma 计算Gamma_us
        # % Gamma_us表示 user u 到 BS s 的SINR
        if math.isnan(pu):
            Pu_temp = self.Pu[user]
        else:
            Pu_temp = pu
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        denominator = 0
        for i in range(0, serverNumber):
            if i != server:
                [Us, n] = self.genUs(x, i)
                if(n > 0):
                    for k in range(0, n):
                        H_temp = self.H
                        if(n == 1):
                            user_temp = Us[0]
                        else:
                            user_temp = Us[k, 0]
                        h_temp = H_temp[user_temp, server, band]
                        denominator = denominator + x[user_temp, i, band] * self.Pu_max * h_temp
        denominator = denominator + self.Sigma_square
        Gamma = Pu_temp * self.H[user, server, band] / denominator
        return Gamma

    def Tcommu(self, x, user, server, band, pu = float('nan')):
        Gamma_us = self.getGamma(x, user, server, band, pu)
        if Gamma_us == 0:
            Tcom = math.inf
        else:
            Tcom = self.Tu_data[user] / self.W / math.log2(1 + Gamma_us)
        return Tcom

    def check_Ttol(self, x):
        [user_vec, server_vec, band_vec] = np.nonzero(x)
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        Pu_max_M = self.Pu_max * np.ones(userNumber)
        Ttol_check_flag = 1
        len_userV = user_vec.size
        for i in range(0, len_userV):
            user = user_vec[i]
            server = server_vec[i]
            band = band_vec[i]
            T_com_Pu_max = self.Tcommu(x, user, server, band, self.Pu_max)
            Ttol_buff = self.Ttol[user, server, band]
            if (T_com_Pu_max > Ttol_buff):
                Ttol_check_flag = 0
                break
        return Ttol_check_flag

    def get_phi(self, lamda_u, beta_time, du, t_local):
        phi_u = lamda_u * beta_time * du / t_local / self.W
        return phi_u

    def get_psi(self, lamda_u, beta_enengy, du, E_local):
        psi_u = lamda_u * beta_enengy * du / E_local / self.W
        return psi_u

    def getupsilong_j(self, x, user, server, band):
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        denominator = 0
        for i in range(0, serverNumber):
            if i != server:
                [Us, n] = self.genUs(x, i)
                for k in range(0, n):
                    if (n == 1):
                        user_temp = Us[0]
                    else:
                        user_temp = Us[k, 0]
                    denominator = denominator + x[user_temp, i, band] * self.Pu_max * self.H[user_temp, server, band]
        denominator = denominator + self.Sigma_square
        upsilong_j = self.H[user, server, band] / denominator
        return upsilong_j

    def getPu_tol(self, du, timeTol, upsilong, sub_bandNumber):
        B = self.W / sub_bandNumber
        buff = 2 ** (du / B / timeTol) - 1
        Pu_tol = buff / upsilong
        return Pu_tol

    def getPi_s(self, x, user, server, band, lamda_u, pu):
        # %GetPi_s 计算Pi_s，论文式(24)
        # % Pi_s表示UPA问题中的计算
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        B = self.W / sub_bandNumber
        du = self.Tu_data[user]
        beta_time_buff = self.beta_time[user]
        beta_enengy_buff = self.beta_enengy[user]

        phi_u = self.get_phi(lamda_u, beta_time_buff, du, self.tu_local[user])
        psi_u = self.get_psi(lamda_u, beta_enengy_buff, du, self.Eu_local[user])

        upsilong_us = self.getupsilong_j(x,user,server,band)
        buff1 = psi_u * math.log2(1 + upsilong_us * pu)
        buff2 = upsilong_us * (phi_u + psi_u * pu) / (1 + upsilong_us * pu) / math.log(2)
        Pi_s = buff1 - buff2
        return Pi_s

    def getGamma_s(self, x, Pu):
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        Gamma_s = 0
        B = self.W / sub_bandNumber
        for server in range(0,serverNumber):
            [Us, n] = self.genUs(x,server)
            if n > 0:
                for user_p in range(0, n):
                    if n == 1:
                        user = Us[0]
                        band = Us[1]
                    else:
                        user = Us[user_p, 0]
                        band = Us[user_p, 1]
                    du = self.Tu_data[user]
                    phi_u = self.get_phi(self.lamda[user],self.beta_time[user], du,self.tu_local[user])
                    psi_u = self.get_psi(self.lamda[user],self.beta_enengy[user],du,self.Eu_local[user])
                    upsilong_us = self.getupsilong_j(x,user,server,band)
                    pu = Pu[user]

                    temp_Gamma_s_u = (phi_u + psi_u * pu)/math.log2(1 + upsilong_us * pu)
                    Gamma_s = Gamma_s + temp_Gamma_s_u
        return Gamma_s

    def UPA_second_order_derivative(self, x):
        # 二分法迭代求Pu(user)，输出为ndarray
        [userNumber, serverNumber, sub_bandNumber] = x.shape
    #     % 是泥土优化，可以二分求解
    # % 首先看在边界处，Pu_min,Pu_max 处取值，如果min处>=0或max处<=0则直接得结果
    # % 如果没得到结果则迭代二分求解，直到迭代结束
        Ttol_M = self.Ttol
        Pu_d = self.Pu_max - self.Pu_min
        Ht = self.H
        Pu = np.zeros(userNumber)
        lambda_vec = self.lamda
        epsilong = Pu_d / 16

        for server in range(0,serverNumber):
            [Us, n] = self.genUs(x,server)
            if n > 0:
                for user_p in range(0, n):
                    if n == 1:
                        user = Us[0]
                        band = Us[1]
                    else:
                        user = Us[user_p, 0]
                        band = Us[user_p, 1]

                    du = self.Tu_data[user]
                    Ttol = Ttol_M[user,server,band]
                    lamda_u = lambda_vec[user]
                    upsilong = self.getupsilong_j(x, user,server,band)
                    Pu_tol = self.getPu_tol(du,Ttol, upsilong, sub_bandNumber)

                    Pi_s_tol = self.getPi_s(x,user,server,band,lamda_u,Pu_tol)
                    Pi_s_max = self.getPi_s(x,user,server,band,lamda_u,self.Pu_max)

                    if(Pi_s_tol >= 0):
                        Pu[user] = Pu_tol
                        continue
                    if(Pi_s_max <= 0):
                        Pu[user] = self.Pu_max
                        continue
                    Pu_up = self.Pu_max
                    Pu_down = Pu_tol
                    while(Pu_up - Pu_down > epsilong):
                        Pu_temp = (Pu_up + Pu_down) / 2
                        Pi_s_temp = self.getPi_s(x,user,server,band,lamda_u,Pu_temp)
                        if(Pi_s_temp <= 0):
                            Pu_down = Pu_temp
                        else:
                            Pu_up = Pu_temp
                    Pu[user] = (Pu_up + Pu_down) / 2
        return Pu

    def RA(self, x):
        # % 根据卸载决策和信道参数计算资源分配决策
        # % 引入网络拓扑时变性，资源分配要满足Ttol
        [F, res_cra] = self.cra(x)  # 计算资源分配
        Jx = 0
        [userNumber, serverNumber, sub_bandNumber] = x.shape
        # % Ttol矩阵表示链路所能维持时长，资源分配决策需要满足Ttol
        # % 要求：当不能满足Ttol的时候，RA输出Jx的值为负数，使得不会选择此方法以提高能效
        check_Ttol_logi = self.check_Ttol(x)
        if (check_Ttol_logi == 0):
            Jx = -1
            Pu = np.zeros(userNumber)
            return Jx, F, Pu

        Pu = self.UPA_second_order_derivative(x)
        Gamma_s = self.getGamma_s(x,Pu)
        for server in range(0, serverNumber):
            [Us, n] = self.genUs(x,server)
            if(n > 0):
                for user in range(0, n):
                    if n == 1:
                        user_p = Us[0]
                    else:
                        user_p = Us[user, 0]
                    Jx = Jx + self.lamda[user_p]

        Jx = Jx - Gamma_s - res_cra
        return Jx, F, Pu



if __name__ == '__main__':
    # print_hi('PyCharm')
    # Op.print_hi('OP!')
    matPath = 'HTtol_check.mat'
    data = io.loadmat(matPath)
    print(data.keys())
    # print(data['H'])
    H = data['H']
    Ttol = data['Ttol']
    # print(type(H))  # H为numpy.ndarray变量
    # print(H.shape)
    [userNumber, serverNumber, sub_bandNumber] = H.shape
    print(H[2, 0, 1])
    Fs = 40e9 * np.ones(serverNumber)  # 服务器运算能力矩阵
    Fu = 1e9 * np.ones(userNumber)  # 用户运算能力矩阵
    # print(Fs)
    Tu_data = 10e5 * np.ones(userNumber)  # 任务数据数
    Tu_circle = 10e9 * np.ones(userNumber)  # 任务运行circle数

    lamda = np.ones(userNumber)  # 用户优先级参数
    beta_time = 0.5 * np.ones(userNumber)  # 优化偏好
    beta_enengy = np.ones(userNumber) - beta_time
    Pu = 0.1 * np.ones(userNumber)
    Pu_max = 1
    Pu_min = 0

    Sigma_square = 1e-18  # 噪声方差
    W = 7000e6  # 系统总带宽7000MHz
    k = 1e-26  # 芯片能耗系数
    # 本地计算时间矩阵
    tu_local = np.zeros(userNumber)
    Eu_local = np.zeros(userNumber)
    for i in range(0, userNumber):
        tu_local[i] = Tu_circle[i] / Fu[i]  # 本地计算时间矩阵
        Eu_local[i] = k * (Fu[i]) ** 2 * Tu_circle[i]  # 本地计算能耗矩阵
    Eta_user = np.zeros(userNumber)
    for i in range(0, userNumber):
        Eta_user[i] = beta_time[i] * Tu_circle[i] * lamda[i] / tu_local[i]

    system_para = paraClass0(beta_time, beta_enengy, Tu_data, Tu_circle, tu_local, Eu_local,
                                  W, H, lamda, Pu, Sigma_square, Fs, Eta_user, Pu_max, Pu_min, Ttol)
    opHJTORA = optimizeAlgorithm(system_para)