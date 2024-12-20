%去除x中的y
function elim=eliminate(x,y)
for n=1:length(y)
    x=x(find(x~=y(n))); %去除x中的y
end
elim=x;