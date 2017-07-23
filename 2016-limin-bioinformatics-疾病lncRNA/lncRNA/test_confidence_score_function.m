clear;
clc;
%159 diseases
disease_GIPSim_matrix=load('disease_GIPSim_matrix.txt');
disease_DOSESim_matrix=load('diseasesim_DOSEwang_matrix.txt');
disease_funSim_matrix=load('diseasesim_funSim_matrix.txt');
disease_icodhprdSim_matrix=load('diseasesim_icodhprd.txt');
disease_jaccardSim_matrix=load('diseasesim_jaccard_go_matrix.txt');
load('DN_Sim');
load('DN_O_Sim');
load('DN_Sim_embedding');
load('DN_O_Sim_embedding');
load('DN_Sim_embedding_nonlabel');
load('DN_O_Sim_embedding_nonlabel');

%117 lncRNAs
lncRNA_GIPSim_matrix=load('lncRNA_GIPSim_matrix.txt');
lncRNA_SeqSim_matrix=load('LncRNA_SeqSim_Matrix.txt');
load('lncN_Sim');
load('lncN_O_Sim');
load('lncN_Sim_embedding');
load('lncN_O_Sim_embedding');

lncRNA_Disease_Matrix=load('lncRNA_Disease_Matrix.txt');

nomarlized_disease_GIPSim_confidence_score_matrix = confidence_score_function(disease_GIPSim_matrix,0);
nomarlized_disease_DOSESim_confidence_score_matrix = confidence_score_function(disease_DOSESim_matrix,0);
nomarlized_disease_funSim_confidence_score_matrix = confidence_score_function(disease_funSim_matrix,0);
nomarlized_disease_icodhprdSim_confidence_score_matrix = confidence_score_function(disease_icodhprdSim_matrix,0);
nomarlized_disease_jaccardSim_confidence_score_matrix = confidence_score_function(disease_jaccardSim_matrix,0);
% [nomarlized_DN_Sim_confidence_score_matrix,AUC_DN] = confidence_score_function(DN_Sim,0);
% [nomarlized_DN_O_Sim_confidence_score_matrix,AUC_DN_O] = confidence_score_function(DN_O_Sim,0);
% nomarlized_DN_Sim_confidence_score_matrix = confidence_score_function(DN_Sim_embedding,0);
% nomarlized_DN_O_Sim_confidence_score_matrix = confidence_score_function(DN_O_Sim_embedding,0);
% nomarlized_DN_Sim_confidence_score_matrix = confidence_score_function(DN_Sim_embedding_nonlabel,0);
% nomarlized_DN_O_Sim_confidence_score_matrix = confidence_score_function(DN_O_Sim_embedding_nonlabel,0);
nomarlized_DN_Sim_confidence_score_matrix = confidence_score_function(DN_Sim_embedding,0);
nomarlized_DN_O_Sim_confidence_score_matrix = confidence_score_function(DN_O_Sim_embedding,0);

nomarlized_lncRNA_GIPSim_confidence_score_matrix = confidence_score_function(lncRNA_GIPSim_matrix,1);
nomarlized_lncRNA_SeqSim_confidence_score_matrix = confidence_score_function(lncRNA_SeqSim_matrix,1);
% [nomarlized_lncN_Sim_confidence_score_matrix,AUC_lncN] = confidence_score_function(lncN_Sim,1);
% [nomarlized_lncN_O_Sim_confidence_score_matrix,AUC_lncN_O] = confidence_score_function(lncN_O_Sim,1);
% [nomarlized_lncN_Sim_confidence_score_matrix,AUC_lncN] = confidence_score_function(lncN_Sim_embedding,1);
% [nomarlized_lncN_O_Sim_confidence_score_matrix,AUC_lncN_O] = confidence_score_function(lncN_O_Sim_embedding,1);
nomarlized_lncN_Sim_confidence_score_matrix = confidence_score_function(lncN_Sim_embedding,1);
nomarlized_lncN_O_Sim_confidence_score_matrix = confidence_score_function(lncN_O_Sim_embedding,1);

% save nomarlized_DN_Sim_confidence_score_matrix
% save nomarlized_DN_O_Sim_confidence_score_matrix
% save nomarlized_lncN_Sim_confidence_score_matrix
% save nomarlized_lncN_O_Sim_confidence_score_matrix
% save lncRNA_Disease_Matrix
save result
