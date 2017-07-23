clear;
clc;

%留一交叉验证
lncRNA_Disease=load('lncRNA_Disease_Matrix.txt');
[M,N]=size(lncRNA_Disease);

% the numbe of links are 285
% l=0;
% for i=1:M
%     for j=1:N
%         if lncRNA_Disease(i,j)==1
%             l=l+1;
%         end
%     end
% end


%lncDN based on neighbor --disease159

indices = crossvalind('Kfold', 117, 10);%将数据样本随机分割为10部分
for i = 1:10 %循环10次，分别取出第i部分作为测试样本，其余两部分作为训练样本
    test = (indices == i);
    train = ~test;
    trainData = lncRNA_Disease(train, :);
    testData = lncRNA_Disease(test, :);
end

% for k=1:10                                    %交叉验证k=5，5个包轮流作为测试集
%     test = (indices == k);                    %获得test1集元素在数据集中对应的单元编号
%     train = ~test;                            %train集元素的编号为非test1元素的编号
%     trainData=lncRNA_Disease;
%     for i=1:length(train)
%         if train(i)==0
%             trainData(i,:)=0;
%         end
%     end 
%     l=0;
%     for i=1:M
%         for j=1:N
%             if trainData(i,j)~=lncRNA_Disease(i,j)
%                 l=l+1;
%             end
%         end
%     end
%     
% end

% % 10折交叉验证
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