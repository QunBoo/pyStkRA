function [Jx, F, Pu] = RA(x,para)
% ����ж�ؾ��ߺ��ŵ�����������Դ�������
% ������������ʱ���ԣ���Դ����Ҫ����Ttol
[F,res_cra] = cra(x,para.Fs,para.Eta_user);        %������Դ����
Jx = 0;                                            %��Ч
[userNumber,serverNumber,sub_bandNumber] = size(x);         %��ȡx��С
% Ttol �����ʾ��·����ά��ʱ������Դ���������Ҫ����Ttol
% Ҫ�󣺵���������Ttol��ʱ��RA���Jx��ֵΪ������ʹ�ò���ѡ��˷����������Ч
check_Ttol_logi = check_Ttol(x,para);
if(check_Ttol_logi == 0)
    Jx = -1;
    Pu = zeros(userNumber,1);
%     Pu = para.Pu * 10^-10;
    return 
end

beta_time = para.beta_time;
beta_enengy = para.beta_enengy;
Pu = UPA_second_order_derivative(x,beta_time,beta_enengy,para,sub_bandNumber,serverNumber);
Gamma_s = getGamma_s(para,x,Pu,serverNumber,sub_bandNumber);
for server = 1:serverNumber
    [Us,n] = genUs(x,server);
    multiplexingNumber = zeros(sub_bandNumber,1);
    for band = 1:sub_bandNumber
        multiplexingNumber(band) = sum(x(:,server,band));
    end
    if(n > 0)
       for user = 1:n
            %%%%%%%%
%             Pi = getPi(x,Us(user,1),server,Us(user,2),sub_bandNumber,multiplexingNumber(Us(user,2)),para.beta_time,para.beta_enengy,para.tu_local,para.Eu_local,para.Tu,para.Pu,para.Ht,para.Sigma_square,para.W);
            Jx = Jx + para.lamda(Us(user,1)) * 1;
       end
    end
end
 Jx = Jx - Gamma_s - res_cra;

end

function Pu = UPA_second_order_derivative(x,beta_time,beta_enengy,para,sub_bandNumber,serverNumber)
% ���ַ�������Pu(user)�����Ϊ����
[userNumber,~,~] = size(x); 
%     Pi_s = getPi_s(x,user,server,band,sub_bandNumber,multiplexingNumber,beta_time,beta_enengy,tu_local,Eu_local,Tu,Pu,Ht,Sigma_square,W,lamda_u,Pu_max,pu);

% �������Ż������Զ������
% ���ȿ��ڱ߽紦��Pu_min,Pu_max ��ȡֵ�����min��>=0��max��<=0��ֱ�ӵý��
% ���û�õ���������������⣬ֱ����������
Ttol_M = para.Ttol;

Pu_d = para.Pu_max - para.Pu_min;
Sigma_square = para.Sigma_square;
Ht = para.Ht;
Pu_max = para.Pu_max;
Pu = zeros(userNumber,1);
tu_local = para.tu_local;
Eu_local = para.Eu_local;
Tu = para.Tu;
lamda_M = para.lamda;
W = para.W;
epsilong = Pu_d / 16;
for server = 1:serverNumber
    [Us,n] = genUs(x,server);
    if(n > 0)
       for user_p = 1:n
            %%%%%%%%
            user = Us(user_p,1);
            band = Us(user_p,2);
            du = Tu(user).data;
            Ttol = Ttol_M(user,server,band);
            lamda_u = lamda_M(user);
            upsilong = getupsilong_j(x,Sigma_square,Ht,user,server,band, Pu_max);
            Pu_tol = getPu_tol(du, W, Ttol, upsilong, sub_bandNumber);
            
            Pi_s_tol = getPi_s(x,user,server,band,sub_bandNumber,1,beta_time,beta_enengy,tu_local,Eu_local,Tu,Pu,Ht,Sigma_square,W,lamda_u,Pu_max,Pu_tol);
            Pi_s_max = getPi_s(x,user,server,band,sub_bandNumber,1,beta_time,beta_enengy,tu_local,Eu_local,Tu,Pu,Ht,Sigma_square,W,lamda_u,Pu_max,Pu_max);
            
            if(Pi_s_tol >= 0)
                Pu(user) = Pu_tol;
                continue;
            end
            if(Pi_s_max <= 0)
                Pu(user) = Pu_max;
                continue;
            end
            Pu_up = Pu_max;
            Pu_down = Pu_tol;
            while(Pu_up - Pu_down > epsilong)
                Pu_temp = (Pu_up + Pu_down)/2;
                Pi_s_temp = getPi_s(x,user,server,band,sub_bandNumber,1,beta_time,beta_enengy,tu_local,Eu_local,Tu,Pu,Ht,Sigma_square,W,lamda_u,Pu_max,Pu_temp);
                if(Pi_s_temp <=0)
                    Pu_down = Pu_temp;
                else
                    Pu_up = Pu_temp;
                end
            end
            Pu(user) = (Pu_up + Pu_down)/2;

       end
    end
end

end

function Gamma_s = getGamma_s(para,x,Pu,serverNumber,sub_bandNumber)
% ����ʽ(23)���������Ǹ�����Դ����������ж��Ч��
% ����ж�ؾ���x����Դ��������Pu���������ж�غ���Դ������ߵĽ��
Gamma_s = 0;
W = para.W;
B = W / sub_bandNumber;
lamda_u = para.lamda;
Sigma_square = para.Sigma_square;
Ht = para.Ht;
Tu = para.Tu;
beta_time = para.beta_time;
beta_enengy = para.beta_enengy;
tu_local = para.tu_local;
Eu_local = para.Eu_local;
Pu_max = para.Pu_max;
for server = 1:serverNumber
   [Us,n] = genUs(x,server);
   if(n > 0)
      for user_p = 1:n
          user = Us(user_p,1);
          band = Us(user_p,2);
          du = Tu(user).data;
          phi_u = get_phi(lamda_u(user),beta_time(user),du,tu_local(user),B);
          psi_u = get_psi(lamda_u(user),beta_enengy(user),du,Eu_local(user),B);
          upsilong_us = getupsilong_j(x,Sigma_square,Ht,user,server,band, Pu_max);
          pu = Pu(user);
          
          temp_Gamma_s_u = (phi_u + psi_u * pu)/log2(1 + upsilong_us * pu);
          Gamma_s = Gamma_s + temp_Gamma_s_u;
      end
   end
end



end
function Pi_s = getPi_s(x,user,server,band,sub_bandNumber,multiplexingNumber,beta_time,beta_enengy,tu_local,Eu_local,Tu,Pu,Ht,Sigma_square,W,lamda_u,Pu_max,pu)
%GetPi_s ����Pi_s������ʽ(24)
% Pi_s��ʾUPA�����еļ���
    B = W / sub_bandNumber;
    du = Tu(user).data;
    beta_time_buff = beta_time(user);
    beta_enengy_buff = beta_enengy(user);
%     lamda_u_i = lamda_u(user);
    phi_u = get_phi(lamda_u,beta_time_buff,du,tu_local(user),B);
    psi_u = get_psi(lamda_u,beta_enengy_buff,du,Eu_local(user),B);
%     upsilong_us ӦΪ�������������
    upsilong_us = getupsilong_j(x,Sigma_square,Ht,user,server,band, Pu_max);
    buff_1 = psi_u * log2(1 + upsilong_us * pu);
    buff_2 = upsilong_us * (phi_u + psi_u * pu) / (1 + upsilong_us * pu) / log(2);
    
    Pi_s = (buff_1 - buff_2);
end
function Pu_tol = getPu_tol(du, W, Ttol, upsilong, sub_bandNumber)
% ���㵥���Ӵ���Pu_tol�����������Ϊ����
    B = W / sub_bandNumber;
    buff = 2^(du/B/Ttol) - 1;
    Pu_tol = buff / upsilong;
end
function upsilong_j = getupsilong_j(G,Sigma_square,H,user,server,band,Pu_max)
    [~,serverNumber,~] = size(G);
    denominator = 0;
    for i = 1:serverNumber
        if i ~= server
            [Us,n] = genUs(G,i);
            for k = 1:n
                denominator = denominator + G(Us(k,1),i,band) * Pu_max * H(Us(k,1),server,band);
            end
        end
    end
    denominator = denominator + Sigma_square;
    upsilong_j = H(user,server,band)/denominator;
end

function phi_u = get_phi(lamda_u,beta_time,du,t_local,W)
    check1 = lamda_u * beta_time ;
    check2 = check1* du/ t_local / W;
    phi_u = lamda_u * beta_time * du / t_local / W;
end

function psi_u = get_psi(lamda_u,beta_enengy,du,E_local,W)
    psi_u = lamda_u * beta_enengy * du / E_local / W;
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
function check_Ttol_logi = check_Ttol(x,para)
    [user_vec,server_vec,band_vec] = ind2sub(size(x),find(x == 1));
    [userNumber,serverNumber,sub_bandNumber] = size(x);         %��ȡx��С
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