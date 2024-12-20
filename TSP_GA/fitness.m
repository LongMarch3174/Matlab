%适应度函数
function f=fitness(fmin,fmax,froad)
f=1-(froad-fmin)/(fmax-fmin);%适应度函数
