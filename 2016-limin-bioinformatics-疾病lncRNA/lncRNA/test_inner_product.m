clear;
clc;

load('DlncN117_embedding.mat');
DlncN117=embedding;
clear embedding
load('lncDN159_embedding.mat');
lncDN159=embedding;
clear embedding
%lncRNAs_diseases
lncRNA_Disease=load('lncRNA_Disease_Matrix.txt');

%inner product
lncN=length(DlncN117);
DN=length(lncDN159);
lncN_Sim=eye(lncN);
DN_Sim=eye(DN);
lncN_DN=zeros(lncN,DN);

for i=1:lncN-1
    for j=i+1:lncN
        lncN_Sim(i,j)=dot(DlncN117(i,:),DlncN117(j,:))/(sqrt(dot(DlncN117(i,:),DlncN117(i,:)))*sqrt(dot(DlncN117(j,:),DlncN117(j,:))));
        lncN_Sim(j,i)=dot(DlncN117(i,:),DlncN117(j,:))/(sqrt(dot(DlncN117(i,:),DlncN117(i,:)))*sqrt(dot(DlncN117(j,:),DlncN117(j,:))));
    end
end

for i=1:DN-1
    for j=i+1:DN
        DN_Sim(i,j)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
        DN_Sim(j,i)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
    end
end

for i=1:lncN
    for j=1:DN
        lncN_DN(i,j)=dot(DlncN117(i,:),lncDN159(j,:));  %/(sqrt(dot(DlncN117(i,:),DlncN117(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))))
    end
end
label=reshape(lncRNA_Disease,[1,18603]);
score=reshape(lncN_DN,[1,18603]);
AUC = roc_curve(label,score,0); 

