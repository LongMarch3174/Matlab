%模拟退火算法求解TSP
clc
clear
close all



% load citys_location.mat  %城市坐标
citys = [81 71;91 3;13 28;91 5;63 10;10 82;28 69;55 32;96 95;96 3;16 44;97 38;96 77;49 80;80 19;
    14 49;42 45;92 65;79 71;96 75;66 28;4 68;85 66;93 16;68 12;76 50;74 96;39 34;66 59;17 22];

x = citys(:,1);%x坐标
y = citys(:,2);%y坐标
m = length(x);%城市个数

d = zeros(m);%距离矩阵存储每个城市之间的距离
for i=1:m-1
    for j=i+1:m
        d(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end
d = d+d';

%获取一个较好初始解
S0=[];
Sum = inf;
for j=1:1000
    S = [1 1+randperm(m-1) 1]; 
    tmp = 0;
    for i = 1:m-1
      tmp = tmp + d(S(i),S(i+1));
    end
 if tmp < Sum
     S0 = S;
     Sum = tmp;
 end
end

%计算初始解的距离
BEST_SUM = fun(d,S);%由于TSP问题需要最后回到起点
Best_S = S;%存储最优路径

%%模拟退火参数
T_begin = 100;%初始温度
T_end = 1;%终止温度
K = 0.99;%温度衰减系数
L = 200;%温度T时的迭代次数

T = T_begin;%开始时迭代温度等于初始温度
while T>T_end
    for i=1:L
        %产生新解
        c = 2+floor((m-1)*rand(1,2));%随机选择两个位置进行交换
        c = sort(c);
        c1 = c(1);c2 = c(2);
        
        new_S = [S(1:c1-1),S(c2:-1:c1),S(c2+1:end)];%新路径
        
        %计算代价函数，即路径差
        df = fun(d,new_S)-fun(d,S);
        
        if df<0
            S = new_S;%接受新解作为下一次迭代的原解
            BEST_SUM = fun(d,new_S);%更新最优距离
            Best_S = S;%更新最优解
        elseif rand<exp(-df/T)
            S = new_S;
        end
    end
    T = K*T;%降温
end

figure(3)
hold on
plot(citys(:,1),citys(:,2),'o','MarkerFaceColor','r')
plot([x(Best_S);x(Best_S(1))],[y(Best_S);y(Best_S(1))])
disp(['最短距离为：',num2str(BEST_SUM)]);
disp('最短路径为：')
disp(Best_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%目标函数（计算路径距离函数）%%%%%%%%%%%%%%%%
function y = fun(d,S)
%d为距离矩阵
%S为路径
%y为根据路径计算得到的距离和

m = length(S);%城市个数

y = 0;
for i=1:m-1
    y = y+d(S(i),S(i+1));
end
end

















