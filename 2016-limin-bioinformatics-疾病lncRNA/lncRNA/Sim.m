clear;
clc;

load('DlncN117_embedding.mat');
DlncN117=embedding;
clear embedding
load('lncDN159-embedding_embedding.mat');
lncDN159=embedding;
clear embedding

lncN=length(DlncN117);
DN=length(lncDN159);
lncN_Sim_embedding=eye(lncN);
DN_Sim_embedding=eye(DN);
lncN_O_Sim_embedding=eye(lncN);
DN_O_Sim_embedding=eye(DN);

for i=1:lncN-1
    for j=i+1:lncN
        lncN_Sim_embedding(i,j)=dot(DlncN117(i,:),DlncN117(j,:))/(sqrt(dot(DlncN117(i,:),DlncN117(i,:)))*sqrt(dot(DlncN117(j,:),DlncN117(j,:))));
        lncN_Sim_embedding(j,i)=dot(DlncN117(i,:),DlncN117(j,:))/(sqrt(dot(DlncN117(i,:),DlncN117(i,:)))*sqrt(dot(DlncN117(j,:),DlncN117(j,:))));
    end
end

for i=1:DN-1
    for j=i+1:DN
        DN_Sim_embedding(i,j)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
        DN_Sim_embedding(j,i)=dot(lncDN159(i,:),lncDN159(j,:))/(sqrt(dot(lncDN159(i,:),lncDN159(i,:)))*sqrt(dot(lncDN159(j,:),lncDN159(j,:))));
    end
end

for i=1:lncN-1
    for j=i+1:lncN
        lncN_O_Sim_embedding(i,j)=norm(DlncN117(i,:)-DlncN117(j,:))/(norm(DlncN117(i,:))*norm(DlncN117(j,:)));
        lncN_O_Sim_embedding(j,i)=norm(DlncN117(i,:)-DlncN117(j,:))/(norm(DlncN117(i,:))*norm(DlncN117(j,:)));
    end
end

for i=1:DN-1
    for j=i+1:DN
        DN_O_Sim_embedding(i,j)=norm(lncDN159(i,:)-lncDN159(j,:))/(norm(lncDN159(i,:))*norm(lncDN159(j,:)));
        DN_O_Sim_embedding(j,i)=norm(lncDN159(i,:)-lncDN159(j,:))/(norm(lncDN159(i,:))*norm(lncDN159(j,:)));
    end
end
save lncN_Sim_embedding
save DN_Sim_embedding
save lncN_O_Sim_embedding
save DN_O_Sim_embedding