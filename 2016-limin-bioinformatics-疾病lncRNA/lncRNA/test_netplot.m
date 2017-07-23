%lncRNAs_diseases
lncRNA_Disease=load('lncRNA_Disease_Matrix.txt');
[m,n]=size(lncRNA_Disease);
A=[eye(m),lncRNA_Disease;lncRNA_Disease',eye(n)];
netplot(A,2);

