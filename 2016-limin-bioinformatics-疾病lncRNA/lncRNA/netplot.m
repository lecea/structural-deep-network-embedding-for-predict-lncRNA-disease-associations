%������netplot
%ʹ�÷���������help netplot
%�޷���ֵ
%����ֻ�ܴ�������ͼ
function netplot(A,flag)
    %���÷�������netplot(A,flag)���޷���ֵ
    %AΪ�ڽӾ�����������
    %flag=1ʱ�����ڽӾ���
    %flag=2ʱ�����������
    %����ֻ�ܴ�������ͼ
    if flag==1      %�ڽӾ����ʾ����ͼ
        ND_netplot(A);
        return;
    end
    
    if flag==2      %���������ʾ����ͼ
        [m n]=size(A);      %����������ڽӾ���
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
