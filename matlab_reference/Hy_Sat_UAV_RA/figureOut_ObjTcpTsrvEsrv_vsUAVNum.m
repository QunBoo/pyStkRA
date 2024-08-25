clear;
close all;

% rng(0);                      %随机数种子
sub_bandNumber = 2;          %子带个数
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';

Max_num = 50;
Min_num = 5;
gap = 5;
UAV_num_x = Min_num:gap:Max_num;
Num_simu_points = (Max_num - Min_num) / gap + 1;
data_show_Obj = zeros(4,Num_simu_points);
data_show_Time = zeros(4,Num_simu_points);
rep_time = 2;

for UAV_num = Min_num:gap:Max_num
    num_begin = UAV_num;
    num_end = UAV_num;
    
    check = stk_aircraft_construct_func(pathStr,num_begin, num_end);
    if(check == 1)
       disp([num2str(UAV_num),' UAVs is constructed'])
    end
% 生成信道增益矩阵以及Ttol矩阵
% 卫星通信信道增益很差，因此考虑任务数据量较小但计算量很大的任务
    [H,Ttol] = stkIriGenGain(sub_bandNumber,pathStr);
    [userNumber, serverNumber, ~ ] = size(H);
    Fs = 40e9 * ones(serverNumber,1);  %服务器运算能力矩阵
Fu = 1e9 * ones(userNumber,1);     %用户运算能力矩阵
T0.data = [];                      %任务数据大小
T0.circle = [];                    %任务所需时钟周期
Tu = repmat(T0,userNumber,1);
Tu_check = 1:userNumber;
for i = 1:userNumber
    if ismember(i,Tu_check)
        Tu(i).data = 10e5;
        Tu(i).circle = 10e9;
    else
        Tu(i).data = 0;
        Tu(i).circle = 0;
    end
end


% 优化偏好，默认为0.5，即能耗偏好和时延偏好相同
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
%用户输出功率矩阵 ———————— UPA问题， 需修改
Pu = 0.0001 * 10^-5 * ones(userNumber,1);
Pu_max = 5;
Pu_min = 0;
Sigma_square = 1e-18;       %噪声方差
W = 7000e6;   %系统总带宽7000MHz
k = 1 * 10^-26;  %芯片能耗系数
%测试算法质量
data_show = cell(5,3);
data_show(1,1) = {'Algorithm'};
data_show(1,2) = {'Objective'};
data_show(1,3) = {'Computing Time'};


% disp('optimize_Exhausted Computing')
% tic;
% [J0, X0, ~] = optimize_exhausted(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
% lamda,Sigma_square,beta_time,beta_enengy,...
% k,...                           % 芯片能耗系数
% userNumber,serverNumber,sub_bandNumber ...
% );
% exhausted_time = toc;
% exhausted_objective = J0;
% data_show(2,1) = {'optimize_Exhausted'};
% data_show(2,2) = {exhausted_objective};
% data_show(2,3) = {exhausted_time};

for i = 1:rep_time
    
% 用户优先级参数，全设定为1
% lamda = rand(userNumber,1)*2;
lamda = ones(userNumber,1);

tic
disp('optimize_greedy Computing')
[J3,X3,F3] = optimize_greedy(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % 芯片能耗系数
userNumber,serverNumber,sub_bandNumber ...
);
greedy_time(i) = toc;
greedy_objective(i) = J3;
data_show(3,1) = {'optimize_Greedy'};
data_show(3,2) = {greedy_objective};
data_show(3,3) = {greedy_time};
data_show_Obj(2,(UAV_num - Min_num)/gap + 1) = mean(greedy_objective);
data_show_Time(2,(UAV_num - Min_num)/gap + 1) = mean(greedy_time);


disp('optimize_stk_hJTORA Computing')
tic;
[J1,X1,F1,Pu1] = optimize_stk_hJTORA(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % 芯片能耗系数
userNumber,serverNumber,sub_bandNumber ...
);
hJTORA_time(i) = toc;
hJTORA_objective(i) = J1;
data_show(4,1) = {'optimize_stk_hJTORA'};
data_show(4,2) = {hJTORA_objective};
data_show(4,3) = {hJTORA_time};
data_show_Obj(3,(UAV_num - Min_num)/gap + 1) = mean(hJTORA_objective);
data_show_Time(3,(UAV_num - Min_num)/gap + 1) = mean(hJTORA_time);

disp('optimize_stk_annealing Computing')
tic;
[J2,X2,F2,Pu2] = optimize_stk_annealing(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % 芯片能耗系数
userNumber,serverNumber,sub_bandNumber,...
10e-9,...                       % 温度下界
0.96,...                        % 温度的下降率
5 ...                           % 邻域解空间的大小
);
annealing_time(i) = toc;
annealing_objective(i) = J2;
data_show(5,1) = {'optimize_stk_annealing'};
data_show(5,2) = {annealing_objective};
data_show(5,3) = {annealing_time};
data_show_Obj(4,(UAV_num - Min_num)/gap + 1) = mean(annealing_objective);
data_show_Time(4,(UAV_num - Min_num)/gap + 1) = mean(annealing_time);

end

[T_temp,E_temp,E_Sum_temp] = time_energyConsumption(X1,F1,Pu1,H,Pu_max,W,Fu,Fs,Tu,...
    lamda,Sigma_square,...
    k);
T_temp_v((UAV_num - Min_num)/gap + 1) = T_temp;
E_temp_v((UAV_num - Min_num)/gap + 1) = E_Sum_temp;

energyConsumElc_v = zeros(size(UAV_num_x));
for i = 1:userNumber    %初始化任务矩阵
        tu_local(i) = Tu(i).circle/Fu(i);   %本地计算时间矩阵
        Eu_local(i) = k * (Fu(i))^2 * Tu(i).circle;    %本地计算能耗矩阵
end
lamda_v = lamda.';
Tu_local_mean = mean(tu_local.*lamda_v);
Eu_local_sum = sum(Eu_local.*lamda_v);
Tu_local_mean_v((UAV_num - Min_num)/gap + 1) = Tu_local_mean;
Eu_local_sum_v((UAV_num - Min_num)/gap + 1)  = Eu_local_sum;
end

figure(1)
hold on
plot(UAV_num_x,data_show_Obj(2,1:(Max_num - Min_num)/gap + 1),'rs-');
plot(UAV_num_x,data_show_Obj(3,1:(Max_num - Min_num)/gap + 1),'mv--');
plot(UAV_num_x,data_show_Obj(4,1:(Max_num - Min_num)/gap + 1),'b*:');
% plot(UAV_num_x,data_show_Obj(1,1:(Max_num - Min_num)/gap + 1),'r^-');

% xlabel('1');
ylabel('目标函数值');
grid on
legend('贪心算法','hJTORA算法','模拟退火算法');
hold off;
% legend('模拟退火算法','hJTORA算法','贪心算法','穷举算法');


figure
hold on
plot(UAV_num_x,data_show_Time(2,1:(Max_num - Min_num)/gap + 1),'rs-');
plot(UAV_num_x,data_show_Time(3,1:(Max_num - Min_num)/gap + 1),'mv--');
plot(UAV_num_x,data_show_Time(4,1:(Max_num - Min_num)/gap + 1),'b*:');
% plot(UAV_num_x,data_show_Time(1,1:(Max_num - Min_num)/gap + 1),'r^-');
set(gca,'yscale','log')
% xlabel('2');
ylabel('计算时间(s)');
grid on
legend('贪心算法','hJTORA算法','模拟退火算法');
% legend('模拟退火算法','hJTORA算法','贪心算法','穷举算法');

figure 
hold on
plot(UAV_num_x,T_temp_v,'mv--');
plot(UAV_num_x,Tu_local_mean_v,'b*:');
legend('启发式算法','本地运算');
ylabel('平均业务时延')
xlabel('用户设备数量');
grid on
hold off


figure 
hold on
plot(UAV_num_x,E_temp_v,'mv--');
% Eu_local_sum_v((UAV_num - Min_num)/gap + 1)  = Eu_local_sum;
plot(UAV_num_x,Eu_local_sum_v,'b*:');
legend('启发式算法','本地运算');
ylabel('系统用户总能耗')
xlabel('用户设备数量');
grid on
hold off
