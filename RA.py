import numpy as np
import paraClass as para
# 废弃代码
def genUs(x,server):
#     生成对应服务器的用户矩阵
    [userNumber,serverNumber,sub_bandNumber] = x.size
    num = 0
    Us = np.array([])
    for user in range(0,userNumber):
        for sub_band in range(0,sub_bandNumber):
            if(x[user,server,sub_band] > 0):
                num = num + 1
                Us_temp = np.array([user, sub_band])
                if(num == 1):
                    Us = Us_temp
                else:
                    Us = np.row_stack((Us, Us_temp))

    return Us, num
def cra(x,Fs,Eta_user):
    # CRA计算资源分配函数
    # x应为np.ndarray变量
    [userNumber,serverNumber,sub_bandNumber] = x.size
    F = np.zeros([userNumber,serverNumber])
    T = 0
    for server in range(0,serverNumber):
        [Us,n] = genUs(x,server)
        if n>0:
            EtaRoot_sum = 0
            for i in range(0,n):
                eta_user_ptemp = Us[i,0]
                eta_user_temp = Eta_user[eta_user_ptemp]
                EtaRoot_sum = EtaRoot_sum + eta_user_temp ** 0.5
            for i in range(0,n):
                eta_user_ptemp = Us[i, 0]
                eta_user_temp = Eta_user[eta_user_ptemp]
                F[eta_user_ptemp,server] = Fs[server] * eta_user_temp ** 0.5 / EtaRoot_sum
            T = T + 1 / Fs(server) * EtaRoot_sum ** 2
    res_cra = T
    return F,res_cra

def RA(x,para0):
    # 根据卸载决策和信道参数计算资源分配决策
    # 引入网络拓扑时变性，资源分配要满足Ttol
    print('Hi,RA')
    [F, res_cra] = cra(x, para0.Fs, para0.Eta_user) #% 计算资源分配
    # print(F)
