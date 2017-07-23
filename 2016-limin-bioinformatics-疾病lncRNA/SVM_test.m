clear;
clc;
%159 diseases
disease_GIPSim_matrix=load('disease_GIPSim_matrix.txt');
disease_DOSESim_matrix=load('diseasesim_DOSEwang_matrix.txt');
disease_funSim_matrix=load('diseasesim_funSim_matrix.txt');
disease_icodhprdSim_matrix=load('diseasesim_icodhprd.txt');
disease_jaccardSim_matrix=load('diseasesim_jaccard_go_matrix.txt');
%117 lncRNAs
lncRNA_GIPSim_matrix=load('lncRNA_GIPSim_matrix.txt');
lncRNA_SeqSim_Matrix=load('LncRNA_SeqSim_Matrix.txt');
%lncRNAs_diseases
lncRNA_Disease_Matrix=load('lncRNA_Disease_Matrix.txt');

[M,N]=size(lncRNA_Disease_Matrix);
disease_confidence_score_matrix=zeros(M,N);
nomarlized_disease_confidence_score_matrix=zeros(M,N);
Max_d=zeros(1,M);
Min_d=1./zeros(1,M);
for i=1:M
    for j=1:N
        numerator=0;
        denominator=0;
        for k=1:N
            if k==j
                continue;
            end
            numerator=numerator+disease_GIPSim_matrix(j,k)*lncRNA_Disease_Matrix(i,k);
            denominator=denominator+disease_GIPSim_matrix(j,k);
        end
        disease_confidence_score_matrix(i,j)=numerator/denominator;
        
        if disease_confidence_score_matrix(i,j)>Max_d(1,i)
            Max_d(1,i)=disease_confidence_score_matrix(i,j);
        end
        if disease_confidence_score_matrix(i,j)<Min_d(1,i)
            Min_d(1,i)=disease_confidence_score_matrix(i,j);
        end
    end
end
for i=1:M
    for j=1:N
        nomarlized_disease_confidence_score_matrix(i,j)=(Max_d(1,i)-disease_confidence_score_matrix(i,j))/(Max_d(1,i)-Min_d(1,i));
    end
end


indices = crossvalind('Kfold', 117, 117);%将数据样本随机分割为117部分
for i = 1:117 %循环117次，分别取出第i部分作为测试样本，其余部分作为训练样本
    test = (indices == i);
    train = ~test;
    trainData = nomarlized_disease_confidence_score_matrix(train, :);
    testData = nomarlized_disease_confidence_score_matrix(test, :);
    trainLabel = lncRNA_Disease_Matrix(train, :);
    testLabel = lncRNA_Disease_Matrix(test, :);
    % SVM网络训练
    model = svmtrain(trainLabel, trainData, '-s 2 -c 1 -g 0.07');
    % SVM网络预测
    [predict_label] = svmpredict(testLabel, testData, model);
    
end
