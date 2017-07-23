clear;
clc;
%159 diseases
disease_GIPSim_matrix=load('disease_GIPSim_matrix.txt');
disease_DOSESim_matrix=load('diseasesim_DOSEwang_matrix.txt');
disease_funSim_matrix=load('diseasesim_funSim_matrix.txt');
disease_icodhprdSim_matrix=load('diseasesim_icodhprd.txt');
disease_jaccardSim_matrix=load('diseasesim_jaccard_go_matrix.txt');
load('lncN_Sim');
load('DN_Sim');
load('lncN_O_Sim');
load('DN_O_Sim');
%117 lncRNAs
lncRNA_GIPSim_matrix=load('lncRNA_GIPSim_matrix.txt');
lncRNA_SeqSim_matrix=load('LncRNA_SeqSim_Matrix.txt');
%lncRNAs_diseases
lncRNA_Disease_Matrix=load('lncRNA_Disease_Matrix.txt');

[M,N]=size(lncRNA_Disease_Matrix);

%disease_GIPSim
disease_GIPSim_confidence_score_matrix=zeros(M,N);
nomarlized_disease_GIPSim_confidence_score_matrix=zeros(M,N);
Max_d_GIPSim=zeros(1,M);
Min_d_GIPSim=1./zeros(1,M);
%disease_DOSESim
disease_DOSESim_confidence_score_matrix=zeros(M,N);
nomarlized_disease_DOSESim_confidence_score_matrix=zeros(M,N);
Max_d_DOSESim=zeros(1,M);
Min_d_DOSESim=1./zeros(1,M);
%disease_funSim
disease_funSim_confidence_score_matrix=zeros(M,N);
nomarlized_disease_funSim_confidence_score_matrix=zeros(M,N);
Max_d_funSim=zeros(1,M);
Min_d_funSim=1./zeros(1,M);
%disease_icodhprdSim
disease_icodhprdSim_confidence_score_matrix=zeros(M,N);
nomarlized_disease_icodhprdSim_confidence_score_matrix=zeros(M,N);
Max_d_icodhprdSim=zeros(1,M);
Min_d_icodhprdSim=1./zeros(1,M);
%disease_jaccardSim
disease_jaccardSim_confidence_score_matrix=zeros(M,N);
nomarlized_disease_jaccardSim_confidence_score_matrix=zeros(M,N);
Max_d_jaccardSim=zeros(1,M);
Min_d_jaccardSim=1./zeros(1,M);

%lncRNA_GIPSim
lncRNA_GIPSim_confidence_score_matrix=zeros(M,N);
nomarlized_lncRNA_GIPSim_confidence_score_matrix=zeros(M,N);
Max_l_GIPSim=zeros(1,N);
Min_l_GIPSim=1./zeros(1,N);

%lncRNA_SeqSim
lncRNA_SeqSim_confidence_score_matrix=zeros(M,N);
nomarlized_lncRNA_SeqSim_confidence_score_matrix=zeros(M,N);
Max_l_SeqSim=zeros(1,N);
Min_l_SeqSim=1./zeros(1,N);

for i=1:M
    for j=1:N
        %disease_GIPSim
        %disease_DOSESim
        %disease_funSim
        %disease_icodhprdSim
        %disease_jaccardSim
        disease_GIPSim_numerator=0;
        disease_GIPSim_denominator=0;
        disease_DOSESim_numerator=0;
        disease_DOSESim_denominator=0;
        disease_funSim_numerator=0;
        disease_funSim_denominator=0;
        disease_icodhprdSim_numerator=0;
        disease_icodhprdSim_denominator=0;
        disease_jaccardSim_numerator=0;
        disease_jaccardSim_denominator=0;
        %lncRNA_GIPSim
        %lncRNA_SeqSim
        lncRNA_GIPSim_numerator=0;
        lncRNA_GIPSim_denominator=0;
        lncRNA_SeqSim_numerator=0;
        lncRNA_SeqSim_denominator=0;
        for k=1:N
            if k==j
                continue;
            end
            %disease_GIPSim
            disease_GIPSim_numerator=disease_GIPSim_numerator+disease_GIPSim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            disease_GIPSim_denominator=disease_GIPSim_denominator+disease_GIPSim_matrix(j,k);
            %disease_DOSESim
            disease_DOSESim_numerator=disease_DOSESim_numerator+disease_DOSESim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            disease_DOSESim_denominator=disease_DOSESim_denominator+disease_DOSESim_matrix(j,k);
            %disease_funSim
            disease_funSim_numerator=disease_funSim_numerator+disease_funSim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            disease_funSim_denominator=disease_funSim_denominator+disease_funSim_matrix(j,k);
            %disease_icodhprdSim
            disease_icodhprdSim_numerator=disease_icodhprdSim_numerator+disease_icodhprdSim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            disease_icodhprdSim_denominator=disease_icodhprdSim_denominator+disease_icodhprdSim_matrix(j,k);
            %disease_jaccardSim
            disease_jaccardSim_numerator=disease_jaccardSim_numerator+disease_jaccardSim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            disease_jaccardSim_denominator=disease_jaccardSim_denominator+disease_jaccardSim_matrix(j,k);
        end
        for k=1:M
            if k==i
                continue;
            end
            %lncRNA_GIPSim
            lncRNA_GIPSim_numerator=lncRNA_GIPSim_numerator+lncRNA_GIPSim_matrix(k,i)*lncRNA_Disease_Matrix(k,j);
            lncRNA_GIPSim_denominator=lncRNA_GIPSim_denominator+lncRNA_GIPSim_matrix(k,i);
            %lncRNA_SeqSim
            lncRNA_SeqSim_numerator=lncRNA_SeqSim_numerator+lncRNA_SeqSim_matrix(k,i)*lncRNA_Disease_Matrix(k,j);
            lncRNA_SeqSim_denominator=lncRNA_SeqSim_denominator+lncRNA_SeqSim_matrix(k,i);
        end
        
        
        %disease_GIPSim
        if disease_GIPSim_denominator==0
            disease_GIPSim_confidence_score_matrix(i,j)=0;
        else
            disease_GIPSim_confidence_score_matrix(i,j)=disease_GIPSim_numerator/disease_GIPSim_denominator;
        end
        if disease_GIPSim_confidence_score_matrix(i,j)>Max_d_GIPSim(1,i)
            Max_d_GIPSim(1,i)=disease_GIPSim_confidence_score_matrix(i,j);
        end
        if disease_GIPSim_confidence_score_matrix(i,j)<Min_d_GIPSim(1,i)
            Min_d_GIPSim(1,i)=disease_GIPSim_confidence_score_matrix(i,j);
        end
        %disease_DOSESim
        if disease_DOSESim_denominator==0
            disease_DOSESim_confidence_score_matrix(i,j)=0;
        else
            disease_DOSESim_confidence_score_matrix(i,j)=disease_DOSESim_numerator/disease_DOSESim_denominator;
        end
        if disease_DOSESim_confidence_score_matrix(i,j)>Max_d_DOSESim(1,i)
            Max_d_DOSESim(1,i)=disease_DOSESim_confidence_score_matrix(i,j);
        end
        if disease_DOSESim_confidence_score_matrix(i,j)<Min_d_DOSESim(1,i)
            Min_d_DOSESim(1,i)=disease_DOSESim_confidence_score_matrix(i,j);
        end
        %disease_funSim
        if disease_funSim_denominator==0
            disease_funSim_confidence_score_matrix(i,j)=0;
        else
            disease_funSim_confidence_score_matrix(i,j)=disease_funSim_numerator/disease_funSim_denominator;
        end
        if disease_funSim_confidence_score_matrix(i,j)>Max_d_funSim(1,i)
            Max_d_funSim(1,i)=disease_funSim_confidence_score_matrix(i,j);
        end
        if disease_funSim_confidence_score_matrix(i,j)<Min_d_funSim(1,i)
            Min_d_funSim(1,i)=disease_funSim_confidence_score_matrix(i,j);
        end
        %disease_icodhprdSim
        if disease_icodhprdSim_denominator==0
            disease_icodhprdSim_confidence_score_matrix(i,j)=0;
        else
            disease_icodhprdSim_confidence_score_matrix(i,j)=disease_icodhprdSim_numerator/disease_icodhprdSim_denominator;
        end
        if disease_icodhprdSim_confidence_score_matrix(i,j)>Max_d_icodhprdSim(1,i)
            Max_d_icodhprdSim(1,i)=disease_icodhprdSim_confidence_score_matrix(i,j);
        end
        if disease_icodhprdSim_confidence_score_matrix(i,j)<Min_d_icodhprdSim(1,i)
            Min_d_icodhprdSim(1,i)=disease_icodhprdSim_confidence_score_matrix(i,j);
        end
        %disease_jaccardSim
        if disease_icodhprdSim_denominator==0
            disease_jaccardSim_confidence_score_matrix(i,j)=0;
        else
            disease_jaccardSim_confidence_score_matrix(i,j)=disease_icodhprdSim_numerator/disease_icodhprdSim_denominator;
        end
        if disease_jaccardSim_confidence_score_matrix(i,j)>Max_d_jaccardSim(1,i)
            Max_d_jaccardSim(1,i)=disease_jaccardSim_confidence_score_matrix(i,j);
        end
        if disease_jaccardSim_confidence_score_matrix(i,j)<Min_d_jaccardSim(1,i)
            Min_d_jaccardSim(1,i)=disease_jaccardSim_confidence_score_matrix(i,j);
        end
        
        %lncRNA_GIPSim
        if lncRNA_GIPSim_denominator==0
            lncRNA_GIPSim_confidence_score_matrix(i,j)=0;
        else
            lncRNA_GIPSim_confidence_score_matrix(i,j)=lncRNA_GIPSim_numerator/lncRNA_GIPSim_denominator;
        end
        if lncRNA_GIPSim_confidence_score_matrix(i,j)>Max_l_GIPSim(1,j)
            Max_l_GIPSim(1,j)=lncRNA_GIPSim_confidence_score_matrix(i,j);
        end
        if lncRNA_GIPSim_confidence_score_matrix(i,j)<Min_l_GIPSim(1,j)
            Min_l_GIPSim(1,j)=lncRNA_GIPSim_confidence_score_matrix(i,j);
        end
        %lncRNA_SeqSim
        if lncRNA_SeqSim_denominator==0
            lncRNA_SeqSim_confidence_score_matrix(i,j)=0;
        else
            lncRNA_SeqSim_confidence_score_matrix(i,j)=lncRNA_SeqSim_numerator/lncRNA_SeqSim_denominator;
        end
        if lncRNA_SeqSim_confidence_score_matrix(i,j)>Max_l_SeqSim(1,j)
            Max_l_SeqSim(1,j)=lncRNA_SeqSim_confidence_score_matrix(i,j);
        end
        if lncRNA_SeqSim_confidence_score_matrix(i,j)<Min_l_SeqSim(1,j)
            Min_l_SeqSim(1,j)=lncRNA_SeqSim_confidence_score_matrix(i,j);
        end
        
    end
end

for i=1:M
    for j=1:N
        if Max_d_GIPSim(1,i)==Min_d_GIPSim(1,i)
            nomarlized_disease_GIPSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_disease_GIPSim_confidence_score_matrix(i,j)=(Max_d_GIPSim(1,i)-disease_GIPSim_confidence_score_matrix(i,j))/(Max_d_GIPSim(1,i)-Min_d_GIPSim(1,i));
        end
        if Max_d_DOSESim(1,i)==Min_d_DOSESim(1,i)
            nomarlized_disease_DOSESim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_disease_DOSESim_confidence_score_matrix(i,j)=(Max_d_DOSESim(1,i)-disease_DOSESim_confidence_score_matrix(i,j))/(Max_d_DOSESim(1,i)-Min_d_DOSESim(1,i));
        end
        if Max_d_funSim(1,i)==Min_d_funSim(1,i)
            nomarlized_disease_funSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_disease_funSim_confidence_score_matrix(i,j)=(Max_d_funSim(1,i)-disease_funSim_confidence_score_matrix(i,j))/(Max_d_funSim(1,i)-Min_d_funSim(1,i));
        end
        if Max_d_icodhprdSim(1,i)==Min_d_icodhprdSim(1,i)
            nomarlized_disease_icodhprdSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_disease_icodhprdSim_confidence_score_matrix(i,j)=(Max_d_icodhprdSim(1,i)-disease_icodhprdSim_confidence_score_matrix(i,j))/(Max_d_icodhprdSim(1,i)-Min_d_icodhprdSim(1,i));
        end
        if Max_d_jaccardSim(1,i)==Min_d_jaccardSim(1,i)
            nomarlized_disease_jaccardSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_disease_jaccardSim_confidence_score_matrix(i,j)=(Max_d_jaccardSim(1,i)-disease_jaccardSim_confidence_score_matrix(i,j))/(Max_d_jaccardSim(1,i)-Min_d_jaccardSim(1,i));  
        end
        if Max_l_GIPSim(1,j)==Min_l_GIPSim(1,j)
            nomarlized_lncRNA_GIPSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_lncRNA_GIPSim_confidence_score_matrix(i,j)=(Max_l_GIPSim(1,j)-lncRNA_GIPSim_confidence_score_matrix(i,j))/(Max_l_GIPSim(1,j)-Min_l_GIPSim(1,j));  
        end
        if Max_l_SeqSim(1,j)==Min_l_SeqSim(1,j)
            nomarlized_lncRNA_SeqSim_confidence_score_matrix(i,j)=0;
        else
            nomarlized_lncRNA_SeqSim_confidence_score_matrix(i,j)=(Max_l_SeqSim(1,j)-lncRNA_SeqSim_confidence_score_matrix(i,j))/(Max_l_SeqSim(1,j)-Min_l_SeqSim(1,j));  
        end
    end
end

% 10折交叉验证
% indices = crossvalind('Kfold', 117, 10);%将数据样本随机分割为10部分
% for i = 1:10 %循环117次，分别取出第i部分作为测试样本，其余部分作为训练样本
%     test = (indices == i);
%     train = ~test;
%     trainData = nomarlized_disease_confidence_score_matrix(train, :);
%     testData = nomarlized_disease_confidence_score_matrix(test, :);
%     trainLabel = lncRNA_Disease_Matrix(train, :);
%     testLabel = lncRNA_Disease_Matrix(test, :);
%     % SVM网络训练
% %     model = svmtrain(trainLabel, trainData, '-s 2 -c 1 -g 0.07');
%     % SVM网络预测
% %     [predict_label] = svmpredict(testLabel, testData, model);
% end

[X_d_GIPSim,Y_d_GIPSim,T_d_GIPSim,AUC_d_GIPSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_disease_GIPSim_confidence_score_matrix,[1,18603]),0);
[X_DOSESim,Y_DOSESim,T_DOSESim,AUC_DOSESim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_disease_DOSESim_confidence_score_matrix,[1,18603]),1);
[X_funSim,Y_funSim,T_funSim,AUC_funSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_disease_funSim_confidence_score_matrix,[1,18603]),1);
[X_icodhprdSim,Y_icodhprdSim,T_icodhprdSim,AUC_icodhprdSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_disease_icodhprdSim_confidence_score_matrix,[1,18603]),1);
[X_jaccardSim,Y_jaccardSim,T_jaccardSim,AUC_jaccardSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_disease_jaccardSim_confidence_score_matrix,[1,18603]),1);
[X_l_GIPSim,Y_l_GIPSim,T_l_GIPSim,AUC_l_GIPSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_lncRNA_GIPSim_confidence_score_matrix,[1,18603]),1);
[X_SeqSim,Y_SeqSim,T_SeqSim,AUC_SeqSim] = perfcurve(reshape(lncRNA_Disease_Matrix,[1,18603]),reshape(nomarlized_lncRNA_SeqSim_confidence_score_matrix,[1,18603]),1);

figure(1);
plot(X_d_GIPSim,Y_d_GIPSim);
hold on;
plot(X_DOSESim,Y_DOSESim);
hold on;
plot(X_funSim,Y_funSim);
hold on;
plot(X_icodhprdSim,Y_icodhprdSim);
hold on;
plot(X_jaccardSim,Y_jaccardSim);
hold on;
plot(X_l_GIPSim,Y_l_GIPSim);
hold on;
plot(X_SeqSim,Y_SeqSim);
hold on;
xlabel('False Positive Rate');
ylabel('True Positive Rate');