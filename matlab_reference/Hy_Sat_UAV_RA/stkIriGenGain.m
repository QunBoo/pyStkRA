function [H, Ttol] = stkIriGenGain(sub_bandNumber,pathStr)
    % 读取scenario
    uiap = actxserver('STK11.application');
    root = uiap.Personality2;
    % UAV场景构建看 stk_aircraft_construct.m
    root.LoadScenario(pathStr);
    sc = root.CurrentScenario;
    % 修改系统时间
    sc.StartTime = '25 Oct 2022 05:50:00.000';
    sc.StopTime = '25 Oct 2022 06:10:00.000';
    check_time = '25 Oct 2022 06:00:00';
    format2datetime = 'dd MMM yyyy HH:mm:ss';
    check_time_date = datetime(check_time,'Locale','en_US','InputFormat',format2datetime);
    % 获取时间内的所有可见性，存储为cell
% 获取所有item的句柄
    allChildren = sc.Children;
    allSatellites = allChildren.GetElements('eSatellite');
    allAircrafts = allChildren.GetElements('eAircraft');
    sat_num = allSatellites.Count;
    uav_num = allAircrafts.Count;
    satNum = sat_num;
    uavNum = uav_num;
    
    sat_name = {};
    uav_name = {};
    
    for i=0:1:sat_num-1
   % eg: sat0=allSatellites.Item(cast(0,'int32'));
        eval(['sat',num2str(i),'=allSatellites.Item(cast(i,''int32''));'])
        eval(['sat_name{i+1,2} = sat',num2str(i),'.instanceName;'])
        eval('sat_name{i+1,1} = i;')
    end
    for i=0:1:uav_num-1
        eval(['uav',num2str(i),'=allAircrafts.Item(cast(i,''int32''));'])
        eval(['uav_name{i+1,2} = uav',num2str(i),'.instanceName;'])
        eval('uav_name{i+1,1} = i;')
    end
    
    % 星机可见矩阵
uav_allAccessIntervals = [];
hWait = waitbar(0,'Please Wait for Computing ASL  ╭(′▽`)╭(′▽`)╯');
for from = 0:1:uav_num-1
   waitbar((from)/uav_num,hWait);
%    from
%        仰角约束条件，由于是加在无人机这个对象上的，所以要加在外层循环，只加一次
   eval(['fromUav=uav',num2str(from),';']);
   uavConstraints = fromUav.AccessConstraints;
   elevation = uavConstraints.AddConstraint('eCstrElevationAngle');
   elevation.EnableMin = 1;
   elevation.Min = 8;
   for to = 0:1:satNum-1
       eval(['toSat=sat',num2str(to),';']);
%        计算可见性
       AS_access = fromUav.GetAccessToObject(toSat);
       AS_access.ComputeAccess();
%        AS_access.get;
       AS_accessIntervals = AS_access.ComputedAccessIntervalTimes;
       if(AS_accessIntervals.Count~=0)
           computedIntervals = AS_accessIntervals.ToArray(0, -1); %时间区间
           temp=cell(AS_accessIntervals.Count,2);                 %存储缓存
           meanRange = cell(AS_accessIntervals.Count,1);
           
           for i=1:1:AS_accessIntervals.Count
               temp{i,1}=from+1;
               temp{i,2}=to+1;
               [str,sto] = AS_accessIntervals.GetInterval(i-1);
               aerDP = AS_access.DataProviders.Item('AER Data').Group.Item('VVLH CBF').Exec(str,sto,1);
               range = cell2mat(aerDP.DataSets.GetDataSetByName('Range').GetValues);
               meanRange{i,1} = mean(range);
           end
           uav_allAccessIntervals = [uav_allAccessIntervals;temp,computedIntervals,meanRange];
           
       end
       
   end
end
close(hWait)
% uav_allAccessIntervals
% 每行表示卫星和无人机之间的连接，第1、2列表示UAV和卫星；第3、4列表示链接维持区间；第5列表示卫星和无人机之间的平均距离
% 第一列是UAV
% 第二列是卫星
[access_row,access_col] = size(uav_allAccessIntervals);
ASL_Access_Cell = cell(access_row, access_col);
for row = 1:access_row
    ASL_Access_Cell{row,1} = uav_allAccessIntervals{row, 1};
    ASL_Access_Cell{row,2} = uav_allAccessIntervals{row, 2};
    buff_cell = uav_allAccessIntervals(row, 3);
    buff_str = buff_cell{1};
    buff_str = buff_str(1:end-4);
    buff_t = datetime(buff_str,'Locale','en_US','InputFormat',format2datetime);
    ASL_Access_Cell{row, 3} = buff_t;
    buff_cell = uav_allAccessIntervals(row, 4);
    buff_str = buff_cell{1};
    buff_str = buff_str(1:end-4);
    buff_t = datetime(buff_str,'Locale','en_US','InputFormat',format2datetime);
    ASL_Access_Cell{row, 4} = buff_t;
    ASL_Access_Cell{row, 5} = uav_allAccessIntervals(row, 5);
end
ASL_range_Ttol = [];
[access_row,access_col] = size(uav_allAccessIntervals);
% 第一列为UAV、第二列为卫星、第三列表示距离，第四列表示剩余时间Ttol，第5列表示信道增益
temp = cell(1,5);
for row = 1:access_row
   if(ASL_Access_Cell{row, 3} < check_time_date)&&(ASL_Access_Cell{row, 4} > check_time_date)
       temp{1} = ASL_Access_Cell{row, 1};
       temp{2} = ASL_Access_Cell{row, 2};
       temp{3} = ASL_Access_Cell{row, 5}{1};
       buff = ASL_Access_Cell{row, 4} - check_time_date;
       temp{4} = seconds(buff);
%        temp{5} = genH(temp{3},6000);             %信道使用C波段6GHz频率
       temp{5} = genH(temp{3},27000);              %信道使用Ka波段17GHz频率
       ASL_range_Ttol = [ASL_range_Ttol; temp];
   end
end
% 将上述计算结果构建为信道矩阵和Ttol矩阵返回
H = zeros(uavNum,satNum,sub_bandNumber);
Ttol = zeros(uavNum,satNum,sub_bandNumber);
[ASL_row,~] = size(ASL_range_Ttol);
for row = 1:ASL_row
    uavIndex = ASL_range_Ttol{row,1};
    satIndex = ASL_range_Ttol{row,2};
    H(uavIndex,satIndex,:) = ASL_range_Ttol{row,5} * ones(1,sub_bandNumber);
    Ttol(uavIndex,satIndex,:) = ASL_range_Ttol{row,4} * ones(1,sub_bandNumber);
end
% 天线增益
H = H * 10^4;
end

function H = genH(range, f)
%     range节点距离(km)； f载波频率（Mhz）
% % 需修改，目前先用地面信道
% gain_DB = 140.7 + 36.7*log10(range/1000);
% c = 0.3 km * MHz

ranges = range;
gain_DB = 20 * log10(f) + 20 * log10(ranges) + 32.4;
H = 1/(10^(gain_DB/10));
end