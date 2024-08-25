function Tcom = Tcommu(x,Pu,H,du,Pu_max,Sigma_square,user,server,band,W)
% ����ж�ؾ��ߺ��ŵ���������ͨ��ʱ��    
% x��ʾж�ؾ��ߣ�Pu��ʾ��ǰ���书�ʣ�H��ʾ�����ŵ����棬du��ʾ��������PumaxΪ����书��
%  Sigma_square�������ʣ�user��server��band��ʾ�û������������ش���W��ʾ����
multiplexingNumber = sum(x(:,server,band));           %���ø���
Gamma_us = getGamma(x,Pu,Sigma_square,H,user,server,band,Pu_max);
Tcom = du / W /log2(1 + Gamma_us) * multiplexingNumber;
end

function Gamma = getGamma(G,Pu,Sigma_square,H,user,server,band,Pu_max)
%GetGamma ����Gamma_us 
% Gamma_us��ʾ user u �� BS s ��SINR
    [~,serverNumber,~] = size(G);
    denominator = 0;
    for i = 1:serverNumber              %���������û������źŵĸ���Ӱ��
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
%GenUs ���ɷ�������Ӧ���û�����
% Us��һ���ʾ�͸÷��������ӵ��û�����2���ʾʹ�õ��ش���
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