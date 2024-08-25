function check = stk_aircraft_construct_func(pathStr,num_begin, num_end)
    uiap = actxserver('STK11.application');
    root = uiap.Personality2;
    % UAV场景构建看 stk_aircraft_construct.m
    root.LoadScenario(pathStr);
    sc = root.CurrentScenario;
    allChildren = sc.Children;
    allAircrafts = allChildren.GetElements('eAircraft');
    uav_num = allAircrafts.Count;
    
    
    
for airc_num = uav_num+1:num_end
airc_num_str = num2str(airc_num);
while(length(airc_num_str)<3)
    airc_num_str = ['0',airc_num_str];
end
aircraft = sc.Children.New('eAircraft',['mycraft',airc_num_str]);

aircraft.SetRouteType('ePropagatorGreatArc');
aircraft.route.Method = 'eDetermineTimeAccFromVel';
%Range_LA：偏离初始点的最大纬度
%Range_LO：偏离初始点的最大经度
Range_LA = 1.0;
Range_LO = 1.0;

%Max_Speed：最大速度，单位：节
%Min_Speed：最小速度，单位：节
Max_Speed = 20;
Min_Speed = 1;

%First_Waypoint：初始位置坐标，纬度、经度
%Number_of_Waypoints：路径点数量
%Starting_Altitude：初始高度，单位：英尺
%Max_Altitude_Change：最大高度变化，单位：英尺

% First_latitute = randi([-89,89]);
% First_longitude = randi([-179,179]);
First_latitute = randi([10,50]);
First_longitude = randi([90,135]);
First_Waypoint = [First_latitute, First_longitude];
Number_of_Waypoints = 60;
Starting_Altitude = 70000;
Max_Altitude_Change = 5000;


Feet_to_Km = 0.3048/1000;
Knots_to_Km_per_Second = 1.852/3600;

% 随机生成1、-1组成的数组，数组长度与路径点数量相等。这个数组用于路径点纬度生成
u = rand(Number_of_Waypoints, 1);
ind = u >= 0.5;
u(ind) = 1;
u(~ind) = -1;

% 随机生成1、-1组成的数组，数组长度与路径点数量相等。这个数组用于路径点经度生成
v = rand(Number_of_Waypoints, 1);
ind = v >= 0.5;
v(ind) = 1;
v(~ind) = -1;

% 随机生成1、-1组成的数组，数组长度与路径点数量相等。这个数组用于路径点高度生成
x = rand(Number_of_Waypoints, 1);
ind = x >= 0.5;
x(ind) = 1;
x(~ind) = -1;

% 随机生成1、-1组成的数组，数组长度与路径点数量相等。这个数组用于路径点速度生成
y = rand(Number_of_Waypoints,1);
ind = y >= 0.5;
y(ind) = 1;
y(~ind) = -1;

%在计算高度时，是按照1000英尺来计算的，所以这里先把最大高度变化量除以1000，后续应用Altitude_Difference变量。
Altitude_Difference = Max_Altitude_Change/1000;

for i=1:Number_of_Waypoints
if i == 1
Results(i,1) = First_Waypoint(1);
Results(i,2) = First_Waypoint(2);
Results(i,3) = (Starting_Altitude)*Feet_to_Km;
Results(i,4) = (randi([Min_Speed,Max_Speed]))*Knots_to_Km_per_Second;
else
Results(i,1) = First_Waypoint(1) + randi(100,1)/100*Range_LA*u(i);
Results(i,2) = First_Waypoint(2) + randi(100,1)/100*Range_LO*v(i);
Results(i,3) = (Starting_Altitude + randi(Altitude_Difference, 1)*1000*x(i))*Feet_to_Km;
Results(i,4) = (randi([Min_Speed,Max_Speed]))*Knots_to_Km_per_Second;
end
end

%添加路径点
for i = 1 : Number_of_Waypoints
    waypoint = aircraft.route.Waypoints.Add();     %句柄
    waypoint.Latitude = Results(i,1);              %纬度
    waypoint.Longitude = Results(i,2);             %经度
    waypoint.Altitude = Results(i,3);              %高度
    waypoint.Speed = Results(i,4);                 %速度
end

aircraft.route.Propagate;
end

root.Save()
check = 1;
end