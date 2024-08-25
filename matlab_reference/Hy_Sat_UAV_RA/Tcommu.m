function Tcom = Tcommu(x,Pu,H,du,Pu_max,Sigma_square,user,server,band,W)
% 根据卸载决策和信道条件计算通信时间    
% x表示卸载决策，Pu表示当前发射功率，H表示整个信道增益，du表示数据量，Pumax为最大发射功率
%  Sigma_square噪声功率，user、server、band表示用户、服务器、载带，W表示带宽
multiplexingNumber = sum(x(:,server,band));           %复用个数
Gamma_us = getGamma(x,Pu,Sigma_square,H,user,server,band,Pu_max);
Tcom = du / W /log2(1 + Gamma_us) * multiplexingNumber;
end

function Gamma = getGamma(G,Pu,Sigma_square,H,user,server,band,Pu_max)
%GetGamma 计算Gamma_us 
% Gamma_us表示 user u 到 BS s 的SINR
    [~,serverNumber,~] = size(G);
    denominator = 0;
    for i = 1:serverNumber              %计算其他用户发射信号的干扰影响
        if i ~= server
            [Us,n] = genUs(G,i);
            for k = 1:n
%                 buff1 = G(Us(k,1),i,band);
%                 buff2 = Pu(Us(k,1));
%                 buff3 = H(Us(k,1),server,band,Pu_max);
                denominator = denominator + G(Us(k,1),i,band) * Pu_max * H(Us(k,1),server,band);
            end
        end
    end
    denominator = denominator + Sigma_square;
    Gamma = Pu(user)*H(user,server,band)/denominator;
end

function [Us,num] = genUs(G,server)
%GenUs 生成服务器对应的用户矩阵
% Us第一项表示和该服务器链接的用户，第2项表示使用的载带。
    [m,~,z] = size(G);
    num = 0;
    Us = [];
    for user = 1:m
        for sub_band = 1:z
            if G(user,server,sub_band) > 0
                num = num + 1;
                Us(num,1) = user;
                Us(num,2) = sub_band;
                break;
            end
        end
    end
end