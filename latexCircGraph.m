function [Y] = latexCircGraph(N,G,order,sX,sY)
n=size(N,1);
angle=2*pi/n;
Y="\scalebox{1}{\begin{tikzpicture}[shorten >=1pt,auto=left,/tikz/initial text=]";
Y=Y+newline+"\def \x {"+sX+"}";
Y=Y+newline+"\def \y {"+sY+"}";
for i=1:n
    if(size(N,2)>1)
        option="[color="+N(order(i),2)+"]";
    end
    Y=Y+newline+"\node ("+num2str(order(i))+")"+option+" at ("+num2str(cos(angle*i))+"*\x,"+num2str(sin(angle*i))+"*\y) {"+N(order(i),1)+"};";
end
Y=Y+newline+"\path[->,every node/.style={font=\scriptsize}]";
for j=1:size(G,1)
    switch G(j,3)
        case 1
            col="blue";
        case 2
            col="red";
        otherwise
            col="black";
    end
    Y=Y+newline+"("+num2str(G(j,1))+") edge["+col+"] ("+num2str(G(j,2))+")";
end
Y=Y+newline+";"+newline+"\end{tikzpicture}}";
clipboard('copy',char(Y))
