clear;
clc;

%��һ������֤
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

indices = crossvalind('Kfold', 117, 10);%��������������ָ�Ϊ10����
for i = 1:10 %ѭ��10�Σ��ֱ�ȡ����i������Ϊ����������������������Ϊѵ������
    test = (indices == i);
    train = ~test;
    trainData = lncRNA_Disease(train, :);
    testData = lncRNA_Disease(test, :);
end

% for k=1:10                                    %������֤k=5��5����������Ϊ���Լ�
%     test = (indices == k);                    %���test1��Ԫ�������ݼ��ж�Ӧ�ĵ�Ԫ���
%     train = ~test;                            %train��Ԫ�صı��Ϊ��test1Ԫ�صı��
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

% % 10�۽�����֤
% indices = crossvalind('Kfold', 117, 10);%��������������ָ�Ϊ10����
% for i = 1:10 %ѭ��117�Σ��ֱ�ȡ����i������Ϊ�������������ಿ����Ϊѵ������
%     test = (indices == i);
%     train = ~test;
%     trainData = nomarlized_disease_confidence_score_matrix(train, :);
%     testData = nomarlized_disease_confidence_score_matrix(test, :);
%     trainLabel = lncRNA_Disease_Matrix(train, :);
%     testLabel = lncRNA_Disease_Matrix(test, :);
%     % SVM����ѵ��
% %     model = svmtrain(trainLabel, trainData, '-s 2 -c 1 -g 0.07');
%     % SVM����Ԥ��
% %     [predict_label] = svmpredict(testLabel, testData, model);
% end