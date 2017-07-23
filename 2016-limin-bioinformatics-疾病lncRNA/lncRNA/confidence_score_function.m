function [nomarlized_confidence_score_matrix,AUC] = confidence_score_function(matrix,flag)
% function nomarlized_confidence_score_matrix = confidence_score_function(matrix,flag)
% caculate nomarlized confidence score of matrix 

%lncRNAs_diseases
lncRNA_Disease_Matrix=load('lncRNA_Disease_Matrix.txt');
[M,N]=size(lncRNA_Disease_Matrix);
confidence_score_matrix=zeros(M,N);
nomarlized_confidence_score_matrix=zeros(M,N);
L=M;
K=N;
if flag==1  %lncRNA
    L=N;
    K=M;
end
Max=zeros(1,L);
Min=1./zeros(1,L);

for i=1:M
    for j=1:N
        numerator=0;
        denominator=0;
        if flag==0
            for k=1:K
                if k==j
                    continue;
                end
                numerator=numerator+matrix(j,k)*lncRNA_Disease_Matrix(i,k);
                denominator=denominator+matrix(j,k);
            end
        else
            for k=1:K
                if k==i
                    continue;
                end
                numerator=numerator+matrix(k,i)*lncRNA_Disease_Matrix(k,j);
                denominator=denominator+matrix(k,i);
            end     
        end
        if denominator==0
            confidence_score_matrix(i,j)=0;
        else
            confidence_score_matrix(i,j)=numerator/denominator;
        end
        if flag==0
            if confidence_score_matrix(i,j)>Max(1,i)
                Max(1,i)=confidence_score_matrix(i,j);
            end
            if confidence_score_matrix(i,j)<Min(1,i)
                Min(1,i)=confidence_score_matrix(i,j);
            end
        else
            if confidence_score_matrix(i,j)>Max(1,j)
                Max(1,j)=confidence_score_matrix(i,j);
            end
            if confidence_score_matrix(i,j)<Min(1,j)
                Min(1,j)=confidence_score_matrix(i,j);
            end
        end   
    end
end

for i=1:M
    for j=1:N
        if flag==0
            if Max(1,i)==Min(1,i)
                nomarlized_confidence_score_matrix(i,j)=0;
            else
                nomarlized_confidence_score_matrix(i,j)=(Max(1,i)-confidence_score_matrix(i,j))/(Max(1,i)-Min(1,i));
            end
        else
            if Max(1,j)==Min(1,j)
                nomarlized_confidence_score_matrix(i,j)=0;
            else
                nomarlized_confidence_score_matrix(i,j)=(Max(1,j)-confidence_score_matrix(i,j))/(Max(1,j)-Min(1,j));  
            end
        end
    end
end

label=reshape(lncRNA_Disease_Matrix,[1,18603]);
score=reshape(nomarlized_confidence_score_matrix,[1,18603]);
AUC = roc_curve(label,score,1); 

% [X,Y,T,AUC] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_confidence_score_matrix,[1,18603]),1);
% disp(['AUC=',num2str(AUC)]);
% figure(1);
% plot(X,Y);
% hold on;
% xlabel('False Positive Rate');
% ylabel('True Positive Rate');
% title('ROC')
