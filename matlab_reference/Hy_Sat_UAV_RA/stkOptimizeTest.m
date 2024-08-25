clear;
close all;

rng(0);                      %���������
sub_bandNumber = 2;          %�Ӵ�����

% �����ŵ���������Լ�Ttol����
% ����ͨ���ŵ�����ܲ��˿���������������С���������ܴ������
pathStr = 'D:\sys\Resource_Allocation\STK\STK\Sc_PGSateNet\PGSateNet.sc';
% num_end = 25;
% num_begin = num_end;
% check = stk_aircraft_construct_func(pathStr,num_begin, num_end);
%     if(check == 1)
%        disp([num2str(num_end),' UAVs is constructed']) 
%     end
% 
% [H,Ttol] = stkIriGenGain(sub_bandNumber,pathStr);

% save HTtol_check_10.mat H Ttol;
load("HTtol_check_15.mat")

[userNumber, serverNumber, ~ ] = size(H);
% userNumber = 10;
% serverNumber = 66;           

Fs = 40e9 * ones(serverNumber,1);  %������������������
Fu = 1e9 * ones(userNumber,1);     %�û�������������

T0.data = [];                      %�������ݴ�С
T0.circle = [];                    %��������ʱ������
Tu = repmat(T0,userNumber,1);
% ��ʼ��������󣬺����ĳ�������������ݺͼ�������
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

% �û����ȼ�������ȫ�趨Ϊ1
lamda = ones(userNumber,1);
% �Ż�ƫ�ã�Ĭ��Ϊ0.5�����ܺ�ƫ�ú�ʱ��ƫ����ͬ
beta_time = 0.5 * ones(userNumber,1);
beta_enengy = ones(userNumber,1) - beta_time;

%�û�������ʾ��� ���������������� UPA���⣬ ���޸�
Pu = 0.1  * ones(userNumber,1);
Pu_max = 1;
Pu_min = 0;

Sigma_square = 1e-18;       %��������
W = 7000e6;   %ϵͳ�ܴ���7000MHz
k = 1 * 10^-26;  %оƬ�ܺ�ϵ��

%�����㷨����
data_show = cell(5,3);
data_show(1,1) = {'Algorithm'};
data_show(1,2) = {'Objective'};
data_show(1,3) = {'Computing Time'};

J0 = 0;
exhausted_time = 0;
disp('optimize_Exhausted Computing')
tic;
% [J0, X0, ~] = optimize_exhausted(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
% lamda,Sigma_square,beta_time,beta_enengy,...
% k,...                           % оƬ�ܺ�ϵ��
% userNumber,serverNumber,sub_bandNumber ...
% );
exhausted_time = toc;
exhausted_objective = J0;
data_show(2,1) = {'optimize_Exhausted'};
data_show(2,2) = {exhausted_objective};
data_show(2,3) = {exhausted_time};
%     

tic
disp('optimize_greedy Computing')
[J3,X3,F3] = optimize_greedy(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % оƬ�ܺ�ϵ��
userNumber,serverNumber,sub_bandNumber ...
);
greedy_time = toc;
greedy_objective = J3;
data_show(3,1) = {'optimize_Greedy'};
data_show(3,2) = {greedy_objective};
data_show(3,3) = {greedy_time};
% 

disp('optimize_stk_hJTORA Computing')
tic;
% J1 = 32.8
[J1,X1,F1,Pu1] = optimize_stk_hJTORA(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
lamda,Sigma_square,beta_time,beta_enengy,...
k,...                           % оƬ�ܺ�ϵ��
userNumber,serverNumber,sub_bandNumber ...
);
hJTORA_time = toc;
hJTORA_objective = J1;
data_show(4,1) = {'optimize_stk_hJTORA'};
data_show(4,2) = {hJTORA_objective};
data_show(4,3) = {hJTORA_time};

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
annealing_time = toc;
annealing_objective = J2;
data_show(5,1) = {'optimize_stk_annealing'};
data_show(5,2) = {annealing_objective};
data_show(5,3) = {annealing_time};

% % ��鸴�����
% check_X3 = sum(X3);
% % reshape(check_X3,serverNumber,sub_bandNumber);
% check_X3_p = find(check_X3 > 1);
% check_X1 = sum(X1);
% % reshape(check_X1,serverNumber,sub_bandNumber);
% check_X1_p = find(check_X1 > 1);
% check_X2 = sum(X2);
% % reshape(check_X2,serverNumber,sub_bandNumber);
% check_X2_p = find(check_X2 > 1);
% 
% % % ��ͼ
figure
% vals_obj = [exhausted_objective,hJTORA_objective,greedy_objective;
%     exhausted_objective,hJTORA_objective,greedy_objective];
vals_obj = [annealing_objective,hJTORA_objective,greedy_objective,exhausted_objective;
    annealing_objective,hJTORA_objective,greedy_objective,exhausted_objective];

b = bar(vals_obj,0.4);
% xlabel('1');
ylabel('Ŀ�꺯��ֵ');
grid on
% legend('����㷨','hJTORA�㷨','̰���㷨');
legend('ģ���˻��㷨','hJTORA�㷨','̰���㷨','����㷨');


figure
% vals_time = [exhausted_time,hJTORA_time,greedy_time;
%     exhausted_time,hJTORA_time,greedy_time];
vals_time = [annealing_time,hJTORA_time,greedy_time,exhausted_time;
    annealing_time,hJTORA_time,greedy_time,exhausted_time];
b = bar(vals_time,0.4);
set(gca,'yscale','log')
% xlabel('2');
ylabel('����ʱ��(s)');
grid on
% legend('����㷨','hJTORA�㷨','̰���㷨');
legend('ģ���˻��㷨','hJTORA�㷨','̰���㷨','����㷨');

[T_temp,E_temp] = time_energyConsumption(X1,F1,Pu1,H,Pu_max,W,Fu,Fs,Tu,...
    lamda,Sigma_square,...
    k)

% disp('Program End');
% vals_obj_4 = vals_obj;
% vals_time_4 = vals_time;
% save vals_save_4_UAVs_Ex.mat vals_obj_4 vals_time_4

