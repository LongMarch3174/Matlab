clc;clear; close all;
% %%%蚁群算法求解TSP问题
% load citys_location.mat  %城市坐标
citys = [81 71;91 3;13 28;91 5;63 10;10 82;28 69;55 32;96 95;96 3;16 44;97 38;96 77;49 80;80 19;
    14 49;42 45;92 65;79 71;96 75;66 28;4 68;85 66;93 16;68 12;76 50;74 96;39 34;66 59;17 22];

x = citys(:,1);%x坐标
y = citys(:,2);%y坐标

n = size(citys,1);%城市个数
%%计算城市之间相互距离
D = zeros(n,n);
for i=1:n-1
   for j=i+1:n
      D(i,j) = sqrt(sum((citys(i,:)-citys(j,:)).^2));
   end
end
D = D+D';
%%初始化参数
m = 6;%蚂蚁个数
alpha = 1;
beta = 5;%启发函数重要程度因子
rho = 0.1;%信息素挥发因子
Q = 10;%常系数
Eta = 1./D;%启发函数
Tau = ones(n,n);%信息素矩阵
Table = zeros(m,n);%路径记录表
iter = 1;%迭代次数初值
iter_max = 200;%最大迭代次数
Route_best = zeros(iter_max,n);%各代最佳路径
Length_best = zeros(iter_max,1);%各代最佳路径的长度
Length_ave = zeros(iter_max,1);%各代路径的平均长度

%%迭代寻找最佳路径
while iter<=iter_max
   %随机产生各个蚂蚁的起点城市
   start = zeros(m,1);
   for i=1:m
       temp = randperm(n);
       start(i) = temp(1);
   end
   Table(:,1) = start;
   %构建解空间
   citys_index = 1:n;
   for i=1:m
       %逐个城市路径选择
       for j=2:n
           tabu = Table(i,1:(j-1));%已访问城市的集合
           allow_index = ~ismember(citys_index,tabu);
           allow = citys_index(allow_index);%待访问城市集合
           P = allow;
           %计算城市转移概率
           for k=1:length(allow)
              P(k) =  Tau(tabu(end),allow(k))^alpha*Eta(tabu(end),allow(k))^beta;
           end
           P = P/sum(P);
           %轮盘赌选择下一个城市
           Pc = cumsum(P);
           target_index = find(Pc>=rand);
           target = allow(target_index(1));
           Table(i,j) = target;
       end
   end
   %计算各个蚂蚁的路径距离
   Length = zeros(m,1);
   for i=1:m
       Route = Table(i,:);
       for j=1:n-1
           Length(i) = Length(i)+D(Route(j),Route(j+1));
       end
       Length(i) = Length(i)+D(Route(n),Route(1));
   end
   %计算最短路径距离及平均距离
   if iter==1
       [min_Length,min_index] = min(Length);
       Length_best(iter) = min_Length;
       Length_ave(iter) = mean(Length);
       Route_best(iter,:) = Table(min_index,:);
   else
       [min_Length,min_index] = min(Length);
       Length_best(iter) = min(Length_best(iter-1),min_Length);
       Length_ave(iter) = mean(Length);
       if Length_best(iter)==min_Length
           Route_best(iter,:) = Table(min_index,:);
       else
           Route_best(iter,:) = Route_best(iter-1,:);
       end
   end
   %更新信息素
   Delta_Tau = zeros(n,n);
   %逐个蚂蚁计算
   for i=1:m
       %逐个城市计算
      for j=1:n-1
          Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1))+Q/Length(i);
      end
      Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1))+Q/Length(i);
   end
   Tau = (1-rho)*Tau+Delta_Tau;
   %迭代次数加1 ，清空路径记录表
   iter = iter+1;
   Table = zeros(m,n);
end

%%结果显示
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['最短距离为：',num2str(Shortest_Length)]);
disp('最短路径为：')
disp([Shortest_Route,Shortest_Route(1)]);

%%绘图
%绘制最优路径图
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
    [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'ro-','MarkerFaceColor','r');

%绘制迭代过程图
figure(2)
plot(Length_ave)
hold on
plot(Length_best)
hold off
legend('每次迭代最优距离','每次迭代平均距离')
xlabel('迭代次数')
ylabel('总路程')
title('蚁群算法迭代收敛图')