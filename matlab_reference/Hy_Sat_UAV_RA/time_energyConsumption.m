function [timeConsum,energyConsum, energySum] = time_energyConsumption(X,F,Pu,H,Pu_max,W,Fu,Fs,Tu,...
    lamda,Sigma_square,...
    k)
[userNumber,serverNumber , subband_Number ] = size(H);
X_sum = zeros(1,userNumber);
p = 1;
% ���ȷ������ж�����
for i = 1:userNumber
   for j = 1:serverNumber
      for ki = 1:subband_Number
          if(X(i,j,ki) > 0)
              X_sum(i) = X_sum(i) + X(i,j,ki); 
              X_temp(p,1) = i;
              X_temp(p,2) = j;
              X_temp(p,3) = ki;
              p = p + 1;
          end
         
      end
   end
end
timeConsum_vec = zeros(1,userNumber);
energyConsum_vec = zeros(1,userNumber);
for i = 1: userNumber
    if(X_sum(i) > 0)
%        ж������µ�ҵ��ʱ��
    row_temp = find(X_temp(:,1) == i);
    user = i;
    server = X_temp(row_temp,2);
    band = X_temp(row_temp,3);
    du = Tu(user).data;
    cu = Tu(user).circle;
    Tcom = Tcommu_Pu(X,Pu,H,du,Pu_max,Sigma_square,user,server,band,W);
    Tcp = Tcomputing(cu,F,user,server,band);
    Etran = Tcom * Pu(user);
    timeConsum_vec(i) = Tcom + Tcp;
    energyConsum_vec(i) = Etran;
    else
%         ���ؼ���ҵ��ʱ��
        Tcp = Tu(i).circle/Fu(i);   %���ؼ���ʱ�����
        Elocal = k * (Fu(i))^2 * Tu(i).circle;    %���ؼ����ܺľ���
        timeConsum_vec(i) = Tcp;
        energyConsum_vec(i) = Elocal;
    end
end

timeConsum = 0;
energyConsum = 0;
for i = 1: userNumber
    timeConsum = timeConsum + lamda(i) * timeConsum_vec(i);
    energyConsum = energyConsum + lamda(i) * energyConsum_vec(i);
end
energySum = energyConsum;
lamda_sum = sum(lamda);
timeConsum = timeConsum / lamda_sum;
energyConsum = energyConsum / lamda_sum;

end

function Tcp = Tcomputing(cu,F,user,server,~)
    F_user_server = F(user,server);
    Tcp = cu / F_user_server;
end
function Tcom = Tcommu_Pu(x,Pu,H,du,Pu_max,Sigma_square,user,server,band,W)
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
                denominator = denominator + G(Us(k,1),i,band) * Pu(Us(k,1)) * H(Us(k,1),server,band);
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