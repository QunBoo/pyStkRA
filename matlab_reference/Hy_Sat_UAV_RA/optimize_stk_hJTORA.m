function [J,X,F,Pu_out] = optimize_stk_hJTORA(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % 芯片能耗系数
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
    para.beta_time = beta_time;               %时间/能耗偏好
    para.beta_enengy = beta_enengy;
    para.Tu = Tu;                             %任务数据
    para.tu_local = tu_local;                 %本地计算时间矩阵
    para.Eu_local = Eu_local;                 %本地计算能耗矩阵
    para.W = W;                               %带宽
    para.Ht = H;                              %信道增益矩阵
    para.lamda = lamda;                       %用户优先级参数
    para.Pu = Pu;                                            %%%%发射功率
    para.Sigma_square = Sigma_square;         %噪声方差
    para.Fs = Fs;                             %服务器运算能力矩阵
    para.Eta_user = Eta_user;                 %CRA所需要的ηu
    para.Pu_max = Pu_max;
    para.Pu_min = Pu_min;
    para.Ttol = Ttol;                        %时变拓扑约束
    
    [J, X, F,Pu_out] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para ...                    % 所需参数
    );
    
end
function [J, X, F,Pu_out] = ta( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    para...                     % 所需参数
)
%TA Task allocation,任务分配算法
%采用论文“Joint Task Offoading and Resource Allocation for Multi-Server Mobile-Edge Computing Networks”的算法

X = genOriginX(userNumber, serverNumber,sub_bandNumber,para);    %得到初始X
% % 注意：初始X卸载要求Ttol>T传输时长————改Fx，在后续UPA问题算法中需要Pu>Ptol
%     根据初始解计算初始目标函数值 J，以及CRA分配F
[J, F,Pu_out] = RA(X,para);

% [T_temp,E_temp] = time_energyConsumption(X,F,Pu_out,para.Ht,para.Pu_max,para.W,1e9 * ones(userNumber,1),para.Fs,para.Tu,...
%     para.lamda,para.Sigma_square,...
%     1 * 10^-26)

    picture = zeros(2,1);
    iterations = 1;
    flag = 1;
while(flag == 1)
        flag = 0;
        [X,J,F,not_find_remove] = remove(X,userNumber,serverNumber,sub_bandNumber,para);
        if not_find_remove == 1
            [X,J,F,not_find_exchange] = exchange(X,userNumber,serverNumber,sub_bandNumber,para);
            if not_find_exchange == 0
                flag = 1;
            end
        end
        picture(iterations,1) = iterations;
        picture(iterations,2) = J;
        iterations = iterations + 1;
end
[J, F,Pu_out] = RA(X,para);
    figure
    plot(picture(:,1),picture(:,2),'b-.');
    title('hJTORA算法进行任务调度优化');
    xlabel('迭代次数');
    ylabel('目标函数值');
end

function [res,old_J,old_F,not_find] = remove(x,userNumber,serverNumber,sub_bandNumber,para)
    user = 1;
    server = 1;
    band = 1;
    not_find = 1;
    [old_J,old_F] = RA(x,para);
    while not_find == 1 && user ~= userNumber && server ~= serverNumber && band ~= sub_bandNumber
        not_find = 1;
        for user=1:userNumber
            for server=1:serverNumber
                for band=1:sub_bandNumber
                    if x(user,server,band) == 1
                        x(user,server,band) = 0;
                        [new_J,new_F] = RA(x,para);
                        if new_J > (1 + 1/1000)*old_J
                            not_find = 0;
                            old_J = new_J;
                            old_F = new_F;
                            break;
                        else
                            x(user,server,band) = 1;
                        end
                    end
                end
                if not_find == 0
                    break
                end
            end
            if not_find == 0
                break
            end
        end
    end
    res = x;
end
function [res,old_J,old_F,not_find] = exchange(x,userNumber,serverNumber,sub_bandNumber,para)
    not_find = 1;
    [old_J,old_F] = RA(x,para);
    x_new = x;
%     [userNumber,serverNumber,sub_bandNumber] = size(x_new);         %提取x大小
    Pu_max_M = para.Pu_max * ones(userNumber,1);
    for user=1:userNumber
        for server=1:serverNumber
            for band=1:sub_bandNumber
                if x(user,server,band) == 0
                    x_new(user,:,:) = 0;
                    x_new(:,server,band) = 0;
                    x_new(user,server,band) = 1;
                    [new_J,new_F] = RA(x_new,para);
                    T_com_Pu_max = Tcommu(x_new,Pu_max_M,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
                    Ttol_check = check_Ttol_s(x_new,para);
                    Ttol_buff = para.Ttol(user,server,band);
                    if (new_J > (1 + 1/1000)*old_J) && (Ttol_check == 1)
%                     if (new_J > (1 + 1/1000)*old_J)
                        not_find = 0;
                        old_J = new_J;
                        old_F = new_F;
                        x = x_new;
                        break;
                    else
                        x_new = x;
                    end
                end
            end
            if not_find == 0
                break
            end
        end
        if not_find == 0
            break
        end
    end
    res = x;
end
function seed = genOriginX(userNumber, serverNumber,sub_bandNumber,para)
%GenLargestSeed
    seed = zeros(userNumber, serverNumber,sub_bandNumber);
    old_J = zeros(userNumber, serverNumber,sub_bandNumber);
    Pu_max_M = para.Pu_max * ones(userNumber,1);
    for user=1:userNumber
        for server=1:serverNumber
            for band=1:sub_bandNumber
                seed(user,server,band) = 1;
                du = para.Tu(user).data;
                W = para.W;
                Ttol_buff = para.Ttol(user,server,band);
                T_com_Pu_max = Tcommu(seed,Pu_max_M,para.Ht,du,para.Pu_max,para.Sigma_square,user,server,band,W);
                if(T_com_Pu_max <= Ttol_buff)
                    
                    [old_J_temp,~] = RA(seed,para);
                    [old_J(user,server,band)] = old_J_temp;
                else
                    old_J(user,server,band) = 0;
                end
                seed(user,server,band) = 0;
            end
        end
    end
    [user,server,band] = ind2sub(size(old_J),find(old_J == max(old_J(:))));
    if(old_J(user(1),server(1),band(1)) > 0)
        seed(user(1),server(1),band(1)) = 1;
    end
end

% function [seed,old_J,F] = genOriginX(userNumber, serverNumber,sub_bandNumber,para)
% X = zeros(userNumber, serverNumber,sub_bandNumber);
% Zeros = zeros(userNumber, serverNumber,sub_bandNumber);
%     for user = 1:userNumber
%        [~,server,~] = ind2sub(size(X),find(para.Ht == max(para.Ht(user,:,1))));
%        if ~isempty(server)
%             sub_band = find(~any(X(:,server(1),:)));
%             if ~isempty(sub_band)
%                 X(user,server(1),sub_band(1)) = 1;
%                 check_Ttol_logi = check_Ttol(X,Zeros,para);
%                 if(check_Ttol_logi == 0)
%                     X(user,server(1),sub_band(1)) = 0;
%                 end
%             end
%        end
%     end
%     [J, F] = RA(X,para);
%     seed = X;
%     old_J = J;
% end


function check_Ttol_logi = check_Ttol(x_new,x_old,para)
    x_diff = x_new - x_old;                               %找出新增的链路
    [user_vec,server_vec,band_vec] = ind2sub(size(x_diff),find(x_diff == 1));
    Ttol_check_flag = 1;
    [userNumber,serverNumber,sub_bandNumber] = size(x_new);         %提取x大小
    Pu_max_M = para.Pu_max * ones(userNumber,1);
    for i = 1 : length(user_vec)
        user = user_vec(i);
        server = server_vec(i);
        band = band_vec(i);
        T_com_Pu_max = Tcommu(x_new,Pu_max_M,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
        Ttol_buff = para.Ttol(user,server,band);
        if(T_com_Pu_max > Ttol_buff)
            Ttol_check_flag = 0;
            break
        end
    end
    check_Ttol_logi = Ttol_check_flag;
end
function check_Ttol_logi = check_Ttol_s(x,para)
    [user_vec,server_vec,band_vec] = ind2sub(size(x),find(x == 1));
    [userNumber,serverNumber,sub_bandNumber] = size(x);         %提取x大小
    Pu_max_M = para.Pu_max * ones(userNumber,1);
    Ttol_check_flag = 1;
    for i = 1 : length(user_vec)
        user = user_vec(i);
        server = server_vec(i);
        band = band_vec(i);
        T_com_Pu_max = Tcommu(x,Pu_max_M,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
        Ttol_buff = para.Ttol(user,server,band);
        if(T_com_Pu_max > Ttol_buff)
            Ttol_check_flag = 0;
            break
        end
    end
    check_Ttol_logi = Ttol_check_flag;
end