function [J, X, F] = optimize_greedy(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
    lamda,Sigma_square,beta_time,beta_enengy,...
    k,...                       % 芯片能耗系数
    userNumber,serverNumber,sub_bandNumber ...
    )

%optimize 负责执行优化操作
    tu_local = zeros(userNumber,1);
    Eu_local = zeros(userNumber,1);
    for i = 1:userNumber    %初始化任务矩阵
        tu_local(i) = Tu(i).circle/Fu(i);   %本地计算时间矩阵
        Eu_local(i) = k * (Fu(i))^2 * Tu(i).circle;    %本地计算能耗矩阵
    end
    Eta_user = zeros(userNumber,1);
    for i=1:userNumber  %计算CRA所需的η
        Eta_user(i) = beta_time(i) * Tu(i).circle * lamda(i) / tu_local(i);
    end
    
    %封装参数
    para.beta_time = beta_time;
    para.beta_enengy = beta_enengy;
    para.Tu = Tu;
    para.tu_local = tu_local;
    para.Eu_local = Eu_local;
    para.W = W;
    para.Ht = H;
    para.lamda = lamda;
    para.Pu = Pu;
    para.Sigma_square = Sigma_square;
    para.Fs = Fs;
    para.Eta_user = Eta_user;
    para.Pu_max = Pu_max;
    para.Pu_min = Pu_min;
    para.Ttol = Ttol;                        %时变拓扑约束
    
   [J, X, F] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para ...                    % 所需参数
    );

end

function [J, X, F] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para ...                    % 所需参数
)
%TA Task allocation,任务分配算法，采用递增判定的贪婪算法

    X = zeros(userNumber, serverNumber,sub_bandNumber);
    for user = 1:userNumber
       [~,server,~] = ind2sub(size(X),find(para.Ht == max(para.Ht(user,:,1))));
       if ~isempty(server)
            sub_band = find(~any(X(:,server(1),:)));
            if ~isempty(sub_band)
                X(user,server(1),sub_band(1)) = 1;
                check_Ttol_logi = check_Ttol(X,para);
                if(check_Ttol_logi == 0)
                    X(user,server(1),sub_band(1)) = 0;
                end
            end
       end
    end
    [J, F] = RA(X,para);
end

function check_Ttol_logi = check_Ttol(x,para)
    [user_vec,server_vec,band_vec] = ind2sub(size(x),find(x == 1));
    Ttol_check_flag = 1;
    for i = 1 : length(user_vec)
        user = user_vec(i);
        server = server_vec(i);
        band = band_vec(i);
        T_com_Pu_max = Tcommu(x,para.Pu,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
        Ttol_buff = para.Ttol(user,server,band);
        if(T_com_Pu_max > Ttol_buff)
            Ttol_check_flag = 0;
            break
        end
    end
    check_Ttol_logi = Ttol_check_flag;
end
