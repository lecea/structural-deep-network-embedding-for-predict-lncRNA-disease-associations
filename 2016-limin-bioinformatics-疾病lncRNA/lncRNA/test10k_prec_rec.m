clear;
clc;

% DN_O=load('nomarlized_DN_O_Sim_confidence_score_matrix.mat');
% DN=load('nomarlized_DN_Sim_confidence_score_matrix.mat');
% lncN_O=load('nomarlized_lncN_O_Sim_confidence_score_matrix.mat');
% lncN=load('nomarlized_lncN_Sim_confidence_score_matrix.mat');
% lncRNA_Disease=load('lncRNA_Disease_Matrix.mat');

Result=load('result.mat');
% print(Result);

DN_O_Sim=Result.nomarlized_DN_O_Sim_confidence_score_matrix;
DN_Sim=Result.nomarlized_DN_Sim_confidence_score_matrix;
lncN_O_Sim=Result.nomarlized_lncN_O_Sim_confidence_score_matrix;
lncN_Sim=Result.nomarlized_lncN_Sim_confidence_score_matrix;
lncRNA_Disease=Result.lncRNA_Disease_Matrix;

disease_DOSESim=Result.nomarlized_disease_DOSESim_confidence_score_matrix;

% label=reshape(lncRNA_Disease,[1,18603]);
% score1=reshape(DN_O_Sim,[1,18603]);
% score2=reshape(DN_Sim,[1,18603]);
% score3=reshape(lncN_O_Sim,[1,18603]);
% score4=reshape(lncN_Sim,[1,18603]);
label=reshape(lncRNA_Disease,[18603,1]);
score1=reshape(DN_O_Sim,[18603,1]);
score2=reshape(DN_Sim,[18603,1]);
score3=reshape(lncN_O_Sim,[18603,1]);
score4=reshape(lncN_Sim,[18603,1]);
score5=reshape(disease_DOSESim,[18603,1]);

%Ωª≤Ê—È÷§



prec_rec(score1, label, 'holdFigure', 1);
% [prec,recc,f1_score, tpr, fpr, thresh]=prec_rec(score1, label, 'holdFigure', 1);
% prec_rec(score2, label, 'holdFigure', 1);
% prec_rec(score3, label, 'holdFigure', 1);
% prec_rec(score4, label, 'holdFigure', 1);
% prec_rec(score5, label, 'holdFigure', 1);
% plot(f1_score);

% figure(2)
% AUC = roc_curve(label,score1,1); 
% hold on;
% AUC = roc_curve(label,score3,1); 
% [stack_x,stack_y,thre,auc]=perfcurve(label_y,deci,a);
% title(['ROC curve of (AUC = ' num2str(auc) ' )']);