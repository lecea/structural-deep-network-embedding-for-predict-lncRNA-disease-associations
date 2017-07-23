% clear;
% clc;
% N=100;
% load('DlncN117.txt');
% % for i = 1:100
% % 	[train,test] = crossvalind('LeaveMOut',N,1);
% % end 
% 
% 
% indices = crossvalind('LeaveMOut',N,1);
% for i = 1:100
%     test = (indices == i); 
%     train = ~test;
% end
clear;
clc;

load('lncDN159-embedding-nonlabel_embedding.mat');
lncDN159=embedding;
clear embedding
DN=length(lncDN159);
DN_Sim_embedding_nonlabel=eye(DN);
DN_O_Sim_embedding_nonlabel=eye(DN);

for i=1:DN-1
    for j=i+1:DN
        DN_Sim_embedding_nonlabel(i,j)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
        DN_Sim_embedding_nonlabel(j,i)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
    end
end

for i=1:DN-1
    for j=i+1:DN
        DN_O_Sim_embedding_nonlabel(i,j)=norm(lncDN159(i,:)-lncDN159(j,:))/(norm(lncDN159(i,:))*norm(lncDN159(j,:)));
        DN_O_Sim_embedding_nonlabel(j,i)=norm(lncDN159(i,:)-lncDN159(j,:))/(norm(lncDN159(i,:))*norm(lncDN159(j,:)));
    end
end
save DN_Sim_embedding_nonlabel
save DN_O_Sim_embedding_nonlabel