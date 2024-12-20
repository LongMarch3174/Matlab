clc;clear; close all;
% %%%��Ⱥ�㷨���TSP����
% load citys_location.mat  %��������
citys = [81 71;91 3;13 28;91 5;63 10;10 82;28 69;55 32;96 95;96 3;16 44;97 38;96 77;49 80;80 19;
    14 49;42 45;92 65;79 71;96 75;66 28;4 68;85 66;93 16;68 12;76 50;74 96;39 34;66 59;17 22];

x = citys(:,1);%x����
y = citys(:,2);%y����

n = size(citys,1);%���и���
%%�������֮���໥����
D = zeros(n,n);
for i=1:n-1
   for j=i+1:n
      D(i,j) = sqrt(sum((citys(i,:)-citys(j,:)).^2));
   end
end
D = D+D';
%%��ʼ������
m = 6;%���ϸ���
alpha = 1;
beta = 5;%����������Ҫ�̶�����
rho = 0.1;%��Ϣ�ػӷ�����
Q = 10;%��ϵ��
Eta = 1./D;%��������
Tau = ones(n,n);%��Ϣ�ؾ���
Table = zeros(m,n);%·����¼��
iter = 1;%����������ֵ
iter_max = 200;%����������
Route_best = zeros(iter_max,n);%�������·��
Length_best = zeros(iter_max,1);%�������·���ĳ���
Length_ave = zeros(iter_max,1);%����·����ƽ������

%%����Ѱ�����·��
while iter<=iter_max
   %��������������ϵ�������
   start = zeros(m,1);
   for i=1:m
       temp = randperm(n);
       start(i) = temp(1);
   end
   Table(:,1) = start;
   %������ռ�
   citys_index = 1:n;
   for i=1:m
       %�������·��ѡ��
       for j=2:n
           tabu = Table(i,1:(j-1));%�ѷ��ʳ��еļ���
           allow_index = ~ismember(citys_index,tabu);
           allow = citys_index(allow_index);%�����ʳ��м���
           P = allow;
           %�������ת�Ƹ���
           for k=1:length(allow)
              P(k) =  Tau(tabu(end),allow(k))^alpha*Eta(tabu(end),allow(k))^beta;
           end
           P = P/sum(P);
           %���̶�ѡ����һ������
           Pc = cumsum(P);
           target_index = find(Pc>=rand);
           target = allow(target_index(1));
           Table(i,j) = target;
       end
   end
   %����������ϵ�·������
   Length = zeros(m,1);
   for i=1:m
       Route = Table(i,:);
       for j=1:n-1
           Length(i) = Length(i)+D(Route(j),Route(j+1));
       end
       Length(i) = Length(i)+D(Route(n),Route(1));
   end
   %�������·�����뼰ƽ������
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
   %������Ϣ��
   Delta_Tau = zeros(n,n);
   %������ϼ���
   for i=1:m
       %������м���
      for j=1:n-1
          Delta_Tau(Table(i,j),Table(i,j+1)) = Delta_Tau(Table(i,j),Table(i,j+1))+Q/Length(i);
      end
      Delta_Tau(Table(i,n),Table(i,1)) = Delta_Tau(Table(i,n),Table(i,1))+Q/Length(i);
   end
   Tau = (1-rho)*Tau+Delta_Tau;
   %����������1 �����·����¼��
   iter = iter+1;
   Table = zeros(m,n);
end

%%�����ʾ
[Shortest_Length,index] = min(Length_best);
Shortest_Route = Route_best(index,:);
disp(['��̾���Ϊ��',num2str(Shortest_Length)]);
disp('���·��Ϊ��')
disp([Shortest_Route,Shortest_Route(1)]);

%%��ͼ
%��������·��ͼ
figure(1)
plot([citys(Shortest_Route,1);citys(Shortest_Route(1),1)],...
    [citys(Shortest_Route,2);citys(Shortest_Route(1),2)],'ro-','MarkerFaceColor','r');

%���Ƶ�������ͼ
figure(2)
plot(Length_ave)
hold on
plot(Length_best)
hold off
legend('ÿ�ε������ž���','ÿ�ε���ƽ������')
xlabel('��������')
ylabel('��·��')
title('��Ⱥ�㷨��������ͼ')