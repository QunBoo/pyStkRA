function [J, X, F] = optimize_exhausted(Fu,Fs,Tu,W,Pu,H,Ttol,Pu_max,Pu_min,...
    lamda,Sigma_square,beta_time,beta_enengy,...
    k,...                       % оƬ�ܺ�ϵ��
    userNumber,serverNumber,sub_bandNumber ...
    )

%optimize ����ִ���Ż�����
    tu_local = zeros(userNumber,1);
    Eu_local = zeros(userNumber,1);
    for i = 1:userNumber    %��ʼ���������
        tu_local(i) = Tu(i).circle/Fu(i);   %���ؼ���ʱ�����
        Eu_local(i) = k * (Fu(i))^2 * Tu(i).circle;    %���ؼ����ܺľ���
    end
    Eta_user = zeros(userNumber,1);
    for i=1:userNumber  %����CRA����Ħ�
        Eta_user(i) = beta_time(i) * Tu(i).circle * lamda(i) / tu_local(i);
    end
    
    %��װ����
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
    para.Ttol = Ttol;                        %ʱ������Լ��
    
   [J, X, F] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para ...                    % �������
    );

end

function [J, X, F] = ta( ...
    userNumber,...              % �û�����
    serverNumber,...            % ����������
    sub_bandNumber,...          % �Ӵ�����
    para ...                    % �������
)
%TA Task allocation,��������㷨��������ٷ�

    x = zeros(userNumber, serverNumber,sub_bandNumber);
    
    global array index;
    
    array = struct;
    index = 1;
    
    search(1,x,userNumber,serverNumber,sub_bandNumber,para);
    
    
    [J,num] = max([array.J]);
    X = array(num).x;
    F = array(num).F;
    
    clear global;
end
 
function search(user,x,userNumber,serverNumber,sub_bandNumber,para)
    global array index;
    Ht = para.Ht;
    Ttol = para.Ttol;
    if user <= userNumber
        for server = 0:serverNumber
            if(user == 1)
                disp("server = " + server + "/33")
            end
            for band = 1:sub_bandNumber
                if(server ~= 0)
                    x(user,server,band) = 1;
                    if(Ht(user,server,band) == 0)
                        x(user,server,band) = 0;
                        continue; 
                    end 
                end
                
                check_Ttol_logi = check_Ttol(x,para);
                if(check_Ttol_logi == 1)
                   [J, F, P] = RA(x,para);
                   array(index).J = J;
                   array(index).F = F;
                   array(index).x = x;
                   index = index + 1; 
                end
                search(user+1,x,userNumber,serverNumber,sub_bandNumber,para);
                if(server ~= 0)
                    x(user,server,band) = 0;
                end
            end
        end
    end
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
