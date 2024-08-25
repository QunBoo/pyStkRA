function check = stk_aircraft_construct_func(pathStr,num_begin, num_end)
    uiap = actxserver('STK11.application');
    root = uiap.Personality2;
    % UAV���������� stk_aircraft_construct.m
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
%Range_LA��ƫ���ʼ������γ��
%Range_LO��ƫ���ʼ�����󾭶�
Range_LA = 1.0;
Range_LO = 1.0;

%Max_Speed������ٶȣ���λ����
%Min_Speed����С�ٶȣ���λ����
Max_Speed = 20;
Min_Speed = 1;

%First_Waypoint����ʼλ�����꣬γ�ȡ�����
%Number_of_Waypoints��·��������
%Starting_Altitude����ʼ�߶ȣ���λ��Ӣ��
%Max_Altitude_Change�����߶ȱ仯����λ��Ӣ��

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

% �������1��-1��ɵ����飬���鳤����·����������ȡ������������·����γ������
u = rand(Number_of_Waypoints, 1);
ind = u >= 0.5;
u(ind) = 1;
u(~ind) = -1;

% �������1��-1��ɵ����飬���鳤����·����������ȡ������������·���㾭������
v = rand(Number_of_Waypoints, 1);
ind = v >= 0.5;
v(ind) = 1;
v(~ind) = -1;

% �������1��-1��ɵ����飬���鳤����·����������ȡ������������·����߶�����
x = rand(Number_of_Waypoints, 1);
ind = x >= 0.5;
x(ind) = 1;
x(~ind) = -1;

% �������1��-1��ɵ����飬���鳤����·����������ȡ������������·�����ٶ�����
y = rand(Number_of_Waypoints,1);
ind = y >= 0.5;
y(ind) = 1;
y(~ind) = -1;

%�ڼ���߶�ʱ���ǰ���1000Ӣ��������ģ����������Ȱ����߶ȱ仯������1000������Ӧ��Altitude_Difference������
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

%���·����
for i = 1 : Number_of_Waypoints
    waypoint = aircraft.route.Waypoints.Add();     %���
    waypoint.Latitude = Results(i,1);              %γ��
    waypoint.Longitude = Results(i,2);             %����
    waypoint.Altitude = Results(i,3);              %�߶�
    waypoint.Speed = Results(i,4);                 %�ٶ�
end

aircraft.route.Propagate;
end

root.Save()
check = 1;
end