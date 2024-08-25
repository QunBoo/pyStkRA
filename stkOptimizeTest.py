
import scipy.io as io
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import paraClass as para
import optimize_stk_hJTORA as hJTORA
import optimize_stk_greedy as greedy



def stkOptimize(matPath):
    # Use a breakpoint in the code line below to debug your script.
    data = io.loadmat(matPath)
    print(data.keys())
    # print(data['H'])
    H = data['H']
    Ttol = data['Ttol']
    # print(type(H))  # H为numpy.ndarray变量
    # print(H.shape)
    [userNumber,serverNumber,sub_bandNumber] = H.shape
    # print(H[2,0,1])
    Fs = 40e9 * np.ones(serverNumber)       #服务器运算能力矩阵
    Fu = 1e9 * np.ones(userNumber)          #用户运算能力矩阵
    # print(Fs)
    Tu_data = 10e5 * np.ones(userNumber)    #任务数据数
    Tu_circle = 10e9 * np.ones(userNumber)  #任务运行circle数

    lamda = np.ones(userNumber)             #用户优先级参数
    beta_time = 0.5 * np.ones(userNumber)   #优化偏好
    beta_enengy = np.ones(userNumber) - beta_time
    # print(beta_enengy)
#     用户输出功率矩阵
    Pu = 0.1 * np.ones(userNumber)
    Pu_max = 1
    Pu_min = 0

    Sigma_square = 1e-18                    #噪声方差
    W = 7000e6                              #系统总带宽7000MHz
    k = 1e-26                                #芯片能耗系数
    # 本地计算时间矩阵
    tu_local = np.zeros(userNumber)
    Eu_local = np.zeros(userNumber)
    for i in range(0, userNumber):
        tu_local[i] = Tu_circle[i] / Fu[i]  # 本地计算时间矩阵
        Eu_local[i] = k * (Fu[i]) ** 2 * Tu_circle[i]  # 本地计算能耗矩阵
    Eta_user = np.zeros(userNumber)
    for i in range(0, userNumber):
        Eta_user[i] = beta_time[i] * Tu_circle[i] * lamda[i] / tu_local[i]

    system_para = para.paraClass0(beta_time, beta_enengy, Tu_data, Tu_circle, tu_local, Eu_local,
                                  W, H, lamda, Pu, Sigma_square, Fs, Eta_user, Pu_max, Pu_min, Ttol)


    #     测试不同算法
    print('optimize_stk_greedy Computing')
    greedyTimeStart = time.time()
    opGreedy = greedy.optimizeGreedy(system_para)
    [X_greedy, J_greedy, F_greedy, Pu_greedy] = opGreedy.ta()
    greedyTimeEnd = time.time()
    Tc_greedy = greedyTimeEnd - greedyTimeStart
    print('greedy Objective Value is', J_greedy)
    print('greedy time cost', Tc_greedy, 's')

    print('optimize_stk_hJTORA Computing')
    hJTORATimeStart = time.time()
    opHJTORA = hJTORA.optimizeHJTORA(system_para)
    [X_hJTORA, J_hJTORA, F_hJTORA, Pu_hJTORA] = opHJTORA.ta()
    hJTORATimeEnd = time.time()
    Tc_hJTORA = hJTORATimeEnd - hJTORATimeStart
    print('hJTORA Objective Value is', J_hJTORA)
    print('hJTORA time cost',Tc_hJTORA,'s')

#     fig
    vals_Obj = np.array([[J_greedy, J_hJTORA],
                         [J_greedy, J_hJTORA]])
    vals_Tc = np.array([[Tc_greedy, Tc_hJTORA],
                        [Tc_greedy, Tc_hJTORA]])
    df1 = pd.DataFrame(vals_Obj, columns=['greedy', 'hJTORA'])
    df1.plot(kind='bar',  grid = True)
    plt.show()
    df2 = pd.DataFrame(vals_Tc, columns=['greedy', 'hJTORA'])
    df2.plot(kind='bar', grid=True)
    plt.yscale('log')
    plt.show()


if __name__ == '__main__':
    plt.close('all')
    # print_hi('PyCharm')
    # Op.print_hi('OP!')
    stkOptimize('HTtol_check_30.mat')