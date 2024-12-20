%ģ���˻��㷨���TSP
clc
clear
close all



% load citys_location.mat  %��������
citys = [81 71;91 3;13 28;91 5;63 10;10 82;28 69;55 32;96 95;96 3;16 44;97 38;96 77;49 80;80 19;
    14 49;42 45;92 65;79 71;96 75;66 28;4 68;85 66;93 16;68 12;76 50;74 96;39 34;66 59;17 22];

x = citys(:,1);%x����
y = citys(:,2);%y����
m = length(x);%���и���

d = zeros(m);%�������洢ÿ������֮��ľ���
for i=1:m-1
    for j=i+1:m
        d(i,j) = sqrt((x(i)-x(j))^2+(y(i)-y(j))^2);
    end
end
d = d+d';

%��ȡһ���Ϻó�ʼ��
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

%�����ʼ��ľ���
BEST_SUM = fun(d,S);%����TSP������Ҫ���ص����
Best_S = S;%�洢����·��

%%ģ���˻����
T_begin = 100;%��ʼ�¶�
T_end = 1;%��ֹ�¶�
K = 0.99;%�¶�˥��ϵ��
L = 200;%�¶�Tʱ�ĵ�������

T = T_begin;%��ʼʱ�����¶ȵ��ڳ�ʼ�¶�
while T>T_end
    for i=1:L
        %�����½�
        c = 2+floor((m-1)*rand(1,2));%���ѡ������λ�ý��н���
        c = sort(c);
        c1 = c(1);c2 = c(2);
        
        new_S = [S(1:c1-1),S(c2:-1:c1),S(c2+1:end)];%��·��
        
        %������ۺ�������·����
        df = fun(d,new_S)-fun(d,S);
        
        if df<0
            S = new_S;%�����½���Ϊ��һ�ε�����ԭ��
            BEST_SUM = fun(d,new_S);%�������ž���
            Best_S = S;%�������Ž�
        elseif rand<exp(-df/T)
            S = new_S;
        end
    end
    T = K*T;%����
end

figure(3)
hold on
plot(citys(:,1),citys(:,2),'o','MarkerFaceColor','r')
plot([x(Best_S);x(Best_S(1))],[y(Best_S);y(Best_S(1))])
disp(['��̾���Ϊ��',num2str(BEST_SUM)]);
disp('���·��Ϊ��')
disp(Best_S);

%%%%%%%%%%%%%%%%%%%%%%%%%%%Ŀ�꺯��������·�����뺯����%%%%%%%%%%%%%%%%
function y = fun(d,S)
%dΪ�������
%SΪ·��
%yΪ����·������õ��ľ����

m = length(S);%���и���

y = 0;
for i=1:m-1
    y = y+d(S(i),S(i+1));
end
end

















