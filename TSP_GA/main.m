clc;
clear all;
%主函数
Map=[82,7;91,38;83,46;71,44;64,60;68,58;83,69;87,76;74,78;71,71];
MaxIter=200;%最大迭代次数
Pop_Size=200;%种群规模
pc=0.8;%选择概率
pm=0.1;  %变异概率
[opt,fval,A]=tspga(Map,MaxIter,Pop_Size,pm,pc);
