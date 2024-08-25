function [J, X, F,Pu_out] = optimize_stk_annealing(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
    lamda,Sigma_square,beta_time,beta_enengy,...
    k,...                       % 芯片能耗系数
    userNumber,serverNumber,sub_bandNumber,...
    T_min,...                   % 温度下界
    alpha,...                   % 温度的下降率
    n ...                      % 邻域解空间的大小
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
    
   [J, X, F, Pu_out] = task_allocation( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    T_min,...                   % 温度下界
    alpha,...                   % 温度的下降率
    n, ...                      % 邻域解空间的大小
    para...                    % 所需参数
    );
end
function [max_objective, X, F, Pu_out] = task_allocation( ...
    userNumber,...              % 用户个数
    serverNumber,...            % 服务器个数
    sub_bandNumber,...          % 子带个数
    T_min,...                   % 温度下界
    alpha,...                   % 温度的下降率
    k, ...                      % 邻域解空间的大小
    para...                     % 所需参数
)
%TA Task allocation,任务分配算法，采用模拟退火算法
[x_old,fx_old,F,Pu_out] = genOriginX(userNumber, serverNumber,sub_bandNumber,para);    %得到初始解
    X = x_old;
    max_objective = fx_old;
%     [max_objective, F, Pu_out] = RA(X,para);
    
    Ht = para.Ht;
    picture = zeros(2,1);
    iterations = 1;
    T = userNumber;
%     max_objective = 0;
    while(T>T_min)
        for I=1:k
            x_new = getneighbourhood(x_old,userNumber, serverNumber,sub_bandNumber,Ht);
            [fx_new, F_new, Pu_new] = RA(x_new,para);
            delta = fx_new-fx_old;
%             T_com_Pu_max = Tcommu(x_new,para.Pu,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
%             Ttol_buff = para.Ttol(user,server,band);
% 检查是否符合Ttol要求
            check_Ttol_logi = check_Ttol(x_new,x_old,para);
            if(check_Ttol_logi == 1)
               if (delta>0)
                x_old = x_new;
                fx_old = fx_new;
                if fx_new > max_objective
                    max_objective = fx_new;
                    X = x_new;
                    F = F_new;
                    Pu_out = Pu_new;
                end
                else
                pro=getProbability(delta,T);
                if(pro>rand)
                    x_old=x_new;
                    fx_old = fx_new;
                end
                end 
            end
        end
        picture(iterations,1) = T;
        picture(iterations,2) = fx_old;
        iterations = iterations + 1;
        T=T*alpha;
    end
%     [max_objective, F, Pu_out] = RA(X,para);
end

function res = getneighbourhood(x,userNumber,serverNumber,sub_bandNumber,Ht)
    
    [user_v,server_v,band_v] = ind2sub(size(x),find(x > 0));
    v_num = length(user_v);
    if v_num == 0
%        disp('v_num == 0! warning!')
%        res = x;
        for i = 1:sub_bandNumber
            user_p = unidrnd(userNumber);
            [~,server,~] = ind2sub(size(x),find(Ht == max(Ht(user_p,:,1))));
            if(isempty(server) == 0)
                subband_p = unidrnd(sub_bandNumber);
                x(user_p,server(1),subband_p) = 1;
            end
            chose_check = rand();
            if(chose_check > 0.5) 
               break; 
            end
        end
        
        
        res = x;
       return;
    else
        user_p = unidrnd(v_num);
        user = user_v(user_p);
        server = server_v(user_p);
        band = band_v(user_p);
    end
    
%     chose_check = rand();
%     if(chose_check > 0.95)
%         user_p = unidrnd(userNumber);
%         [~,server_temp,~] = ind2sub(size(x),find(Ht == max(Ht(user_p,:,1))));
%         if(isempty(server_temp) == 0)
% %                 subband_p = unidrnd(sub_bandNumber);
%                 subband_p = find(~any(x(:,server_temp(1),:)));
%                 if ~isempty(subband_p)
%                     band_num = length(subband_p);
%                     chose_p = unidrnd(band_num);
%                 x(user_p,server_temp(1),subband_p(chose_p)) = 1;
%                 end
%         end
%     end
    
    %两种扰动方式，交换或者赋值
    chosen = rand;
    if chosen > 0.2
        if chosen < 0.85   %55%的概率更改用户的服务器（选择offload）
            x(user,server,band) = 0;
            vary_server = unidrnd(serverNumber);    %目标服务器
            check_emp = sum(x(:,vary_server,:));
            check_emp2 = reshape(check_emp,1,sub_bandNumber);
            emp_p = find(check_emp2 == 0);
            if(isempty(emp_p) == 0)
                vary_band=emp_p(randperm(numel(emp_p),1));
                x(user,vary_server,vary_band) = 1;
            else
                x(user,server,band) = 1;
            end

            
        else    %25%的概率更改用户的频带（选择offload）
            if sub_bandNumber ~= 1
                x(user,server,band) = 0;
                check_emp = sum(x(:,server,:));
                check_emp2 = reshape(check_emp,1,sub_bandNumber);
                emp_p = find(check_emp2 == 0);
                if(isempty(emp_p) == 0)
                    vary_band=emp_p(randperm(numel(emp_p),1));
                    x(user,server,vary_band) = 1;
                else
                    x(user,server,band) = 1;
                end
%                 vary_band = emp_p(randperm(numel(emp_p),1));    %目标频带
% %                 while vary_band == band
% %                     vary_band = unidrnd(sub_bandNumber);
% %                 end
%                 x(user,server,vary_band) = 1;
            end
        end
    else 
        if chosen > 0.05  %15%的概率交换两个用户的服务器和频带
            if userNumber ~= 1
                user_other = unidrnd(userNumber);    %指定另一个用户
                while user_other == user
                    user_other = unidrnd(userNumber);
                end
                flag_found = 0;
                for server_other = 1:serverNumber
                    for band_other=1:sub_bandNumber
                        if x(user_other,server_other,band_other) ~= 0
                            flag_found = 1;
                            break;  %找到另一个用户所分配的服务器和频带
                        end
                    end
                    if flag_found == 1
                        break;
                    end
                end
                xValue =  x(user,server,band);
                xValue_other =  x(user_other,server_other,band_other);
                x(user,server,band) = 0;
                x(user_other,server_other,band_other) = 0;
                x(user,server_other,band_other) = xValue_other;  %更改频带和服务器
                x(user_other,server,band) = xValue;
            end
        else    %5%的概率改变该用户的决策
            x(user,server,band) = 1 - x(user,server,band);
        end
    end
    res = x;
% % % % % % % % %     check
%     check_X2 = sum(x);
% % reshape(check_X2,serverNumber,sub_bandNumber);
% check_X2_p = find(check_X2 > 1);
% if(isempty(check_X2_p) == 0)
%    disp('warning!')
%    res;
% end
end

function p = getProbability(delta,t)
    p = exp(delta/t);
end

% function [seed,J,F,Pu] = genOriginX(userNumber, serverNumber,sub_bandNumber,para)
% %GenLargestSeed
%     seed = zeros(userNumber, serverNumber,sub_bandNumber);
%     old_J = zeros(userNumber, serverNumber,sub_bandNumber);
%     for user=1:userNumber
%         for server=1:serverNumber
%             for band=1:sub_bandNumber
%                 seed(user,server,band) = 1;
%                 du = para.Tu(user).data;
%                 W = para.W;
%                 Ttol_buff = para.Ttol(user,server,band);
%                 T_com_Pu_max = Tcommu(seed,para.Pu,para.Ht,du,para.Pu_max,para.Sigma_square,user,server,band,W);
%                 if(T_com_Pu_max <= Ttol_buff)
%                     [old_J(user,server,band),~] = RA(seed,para);
%                 else
%                     old_J(user,server,band) = 0;
%                 end
%                 seed(user,server,band) = 0;
%             end
%         end
%     end
%     [user,server,band] = ind2sub(size(old_J),find(old_J == max(old_J(:))));
%     if(old_J(user(1),server(1),band(1)) > 0)
%         seed(user(1),server(1),band(1)) = 1;
%     end
%     [J, F,Pu] = RA(seed,para);
% end

function [seed,J,F,Pu] = genOriginX(userNumber, serverNumber,sub_bandNumber,para)
%GenRandSeed    生成满足约束的随机种子矩阵
    seed = zeros(userNumber, serverNumber,sub_bandNumber);
    old_J = 0;
    F = zeros(userNumber,serverNumber);
    for user=1:userNumber
        find = 0;
        for server=1:serverNumber
            for band=1:sub_bandNumber
                seed(user,server,band) = 1;
                du = para.Tu(user).data;
                W = para.W;
                Ttol_buff = para.Ttol(user,server,band);
                T_com_Pu_max = Tcommu(seed,para.Pu,para.Ht,du,para.Pu_max,para.Sigma_square,user,server,band,W);
                [new_J,new_F] = RA(seed,para);
                sum_seed = sum(seed);
                check_temp = 0;
                for p_server = 1:serverNumber
                   for p_band = 1:sub_bandNumber
                      if( sum_seed(1,p_server,p_band) > 1)
                          check_temp = 1;
                      end
                   end
                end
                if (new_J > old_J) && (T_com_Pu_max <= Ttol_buff)&&(check_temp == 0)
                    old_J = new_J;
                    F = new_F;
                    find = 1;
                    break;
                else
                    seed(user,server,band) = 0;
                end
            end
            if find == 1
                break;
            end
        end
    end
    [J, F,Pu] = RA(seed,para);
end

function check_Ttol_logi = check_Ttol(x_new,x_old,para)
    x_diff = x_new - x_old;                               %找出新增的链路
    [user_vec,server_vec,band_vec] = ind2sub(size(x_diff),find(x_diff == 1));
    Ttol_check_flag = 1;
    for i = 1 : length(user_vec)
        user = user_vec(i);
        server = server_vec(i);
        band = band_vec(i);
        T_com_Pu_max = Tcommu(x_new,para.Pu,para.Ht,para.Tu(user).data,para.Pu_max,para.Sigma_square,user,server,band,para.W);
        Ttol_buff = para.Ttol(user,server,band);
        if(T_com_Pu_max > Ttol_buff)
            Ttol_check_flag = 0;
            break
        end
    end
    check_Ttol_logi = Ttol_check_flag;
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