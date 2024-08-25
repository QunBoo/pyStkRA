clear;
close all;

% rng(0);                      %���������
sub_bandNumber = 2;          %�Ӵ�����
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
% �����ŵ���������Լ�Ttol����
% ����ͨ���ŵ�����ܲ��˿���������������С���������ܴ������
    [H,Ttol] = stkIriGenGain(sub_bandNumber,pathStr);
    [userNumber, serverNumber, ~ ] = size(H);
    Fs = 40e9 * ones(serverNumber,1);  %������������������
Fu = 1e9 * ones(userNumber,1);     %�û�������������
T0.data = [];                      %�������ݴ�С
T0.circle = [];                    %��������ʱ������
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


% �Ż�ƫ�ã�Ĭ��Ϊ0.5�����ܺ�ƫ�ú�ʱ��ƫ����ͬ
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;
%�û�������ʾ��� ���������������� UPA���⣬ ���޸�
Pu = 0.0001 * 10^-5 * ones(userNumber,1);
Pu_max = 5;
Pu_min = 0;
Sigma_square = 1e-18;       %��������
W = 7000e6;   %ϵͳ�ܴ���7000MHz
k = 1 * 10^-26;  %оƬ�ܺ�ϵ��
%�����㷨����
data_show = cell(5,3);
data_show(1,1) = {'Algorithm'};
data_show(1,2) = {'Objective'};
data_show(1,3) = {'Computing Time'};


% disp('optimize_Exhausted Computing')
% tic;
% [J0, X0, ~] = optimize_exhausted(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
% lamda,Sigma_square,beta_time,beta_enengy,...
% k,...                           % оƬ�ܺ�ϵ��
% userNumber,serverNumber,sub_bandNumber ...
% );
% exhausted_time = toc;
% exhausted_objective = J0;
% data_show(2,1) = {'optimize_Exhausted'};
% data_show(2,2) = {exhausted_objective};
% data_show(2,3) = {exhausted_time};

for i = 1:rep_time
    
% �û����ȼ�������ȫ�趨Ϊ1
% lamda = rand(userNumber,1)*2;
lamda = ones(userNumber,1);

tic
disp('optimize_greedy Computing')
[J3,X3,F3] = optimize_greedy(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % оƬ�ܺ�ϵ��
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
k,...                           % оƬ�ܺ�ϵ��
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
k,...                           % оƬ�ܺ�ϵ��
userNumber,serverNumber,sub_bandNumber,...
10e-9,...                       % �¶��½�
0.96,...                        % �¶ȵ��½���
5 ...                           % �����ռ�Ĵ�С
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
for i = 1:userNumber    %��ʼ���������
        tu_local(i) = Tu(i).circle/Fu(i);   %���ؼ���ʱ�����
        Eu_local(i) = k * (Fu(i))^2 * Tu(i).circle;    %���ؼ����ܺľ���
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
ylabel('Ŀ�꺯��ֵ');
grid on
legend('̰���㷨','hJTORA�㷨','ģ���˻��㷨');
hold off;
% legend('ģ���˻��㷨','hJTORA�㷨','̰���㷨','����㷨');


figure
hold on
plot(UAV_num_x,data_show_Time(2,1:(Max_num - Min_num)/gap + 1),'rs-');
plot(UAV_num_x,data_show_Time(3,1:(Max_num - Min_num)/gap + 1),'mv--');
plot(UAV_num_x,data_show_Time(4,1:(Max_num - Min_num)/gap + 1),'b*:');
% plot(UAV_num_x,data_show_Time(1,1:(Max_num - Min_num)/gap + 1),'r^-');
set(gca,'yscale','log')
% xlabel('2');
ylabel('����ʱ��(s)');
grid on
legend('̰���㷨','hJTORA�㷨','ģ���˻��㷨');
% legend('ģ���˻��㷨','hJTORA�㷨','̰���㷨','����㷨');

figure 
hold on
plot(UAV_num_x,T_temp_v,'mv--');
plot(UAV_num_x,Tu_local_mean_v,'b*:');
legend('����ʽ�㷨','��������');
ylabel('ƽ��ҵ��ʱ��')
xlabel('�û��豸����');
grid on
hold off


figure 
hold on
plot(UAV_num_x,E_temp_v,'mv--');
% Eu_local_sum_v((UAV_num - Min_num)/gap + 1)  = Eu_local_sum;
plot(UAV_num_x,Eu_local_sum_v,'b*:');
legend('����ʽ�㷨','��������');
ylabel('ϵͳ�û����ܺ�')
xlabel('�û��豸����');
grid on
hold off
