function [F,T] = cra(G,Fs,Eta_user)
%CRA computer resourses allocation 运算能力分配
    [userNumber,serverNumber,~] = size(G);
    F = zeros(userNumber,serverNumber);
    T = 0;
    for server = 1:serverNumber
        [Us,n] = genUs(G,server);
        if n > 0 
            EtaRoot_sum = sum(Eta_user(Us(:,1)).^(0.5));
            F(Us(:,1),server) = Fs(server) * Eta_user(Us(:,1)).^(0.5) / EtaRoot_sum;
            T = T + 1/Fs(server) * EtaRoot_sum^2;
        end
    end
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