%函数名netplot
%使用方法输入请help netplot
%无返回值
%函数只能处理无向图
function netplot(A,flag)
    %调用方法输入netplot(A,flag)，无返回值
    %A为邻接矩阵或关联矩阵
    %flag=1时处理邻接矩阵
    %flag=2时处理关联矩阵
    %函数只能处理无向图
    if flag==1      %邻接矩阵表示无向图
        ND_netplot(A);
        return;
    end
    
    if flag==2      %关联矩阵表示无向图
        [m n]=size(A);      %关联矩阵变邻接矩阵
        W=zeros(m,m);
        for i=1:n
            a=find(A(:,i)~=0);
            W(a(1),a(2))=1;
            W(a(2),a(1))=1;
        end
        ND_netplot(W);
        return;
    end 
end
