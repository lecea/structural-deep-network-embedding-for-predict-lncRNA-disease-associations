import random
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab 
import seaborn as sns
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

#set background
#sns.set(style="white",color_codes=True)
num_bins = 50
#load data
disease_GIPSim_matrix=np.loadtxt('disease_GIPSim_matrix.txt')
diseasesim_DOSEwang_matrix=np.loadtxt('diseasesim_DOSEwang_matrix.txt')
diseasesim_funSim_matrix=np.loadtxt('diseasesim_funSim_matrix.txt')
diseasesim_icodhprd=np.loadtxt('diseasesim_icodhprd.txt')
diseasesim_jaccard_go_matrix=np.loadtxt('diseasesim_jaccard_go_matrix.txt')
lncRNA_GIPSim_matrix=np.loadtxt('lncRNA_GIPSim_matrix.txt')
LncRNA_SeqSim_Matrix=np.loadtxt('LncRNA_SeqSim_Matrix.txt')
#disease_GIPSim_matrix=tril(disease_GIPSim_matrix,-1)

disease_GIPSim=disease_GIPSim_matrix.flatten()
diseasesim_DOSEwang=diseasesim_DOSEwang_matrix.flatten()
diseasesim_funSim=diseasesim_funSim_matrix.flatten()
diseasesim_icodhprd=diseasesim_icodhprd.flatten()
diseasesim_jaccard=diseasesim_jaccard_go_matrix.flatten()
lncRNA_GIPSim=lncRNA_GIPSim_matrix.flatten()
LncRNA_SeqSim=LncRNA_SeqSim_Matrix.flatten()

disease_GIPSim_mu=np.mean(disease_GIPSim)#mean
disease_GIPSim_sigma=np.std(disease_GIPSim,ddof=1)
diseasesim_DOSEwang_mu=np.mean(diseasesim_DOSEwang)#mean
diseasesim_DOSEwang_sigma=np.std(diseasesim_DOSEwang,ddof=1)
diseasesim_funSim_mu=np.mean(diseasesim_funSim)#mean
diseasesim_funSim_sigma=np.std(diseasesim_funSim,ddof=1)
diseasesim_icodhprd_mu=np.mean(diseasesim_icodhprd)#mean
diseasesim_icodhprd_sigma=np.std(diseasesim_icodhprd,ddof=1)
diseasesim_jaccard_mu=np.mean(diseasesim_jaccard)#mean
diseasesim_jaccard_sigma=np.std(diseasesim_jaccard,ddof=1)
lncRNA_GIPSim_mu=np.mean(lncRNA_GIPSim)#mean
lncRNA_GIPSim_sigma=np.std(lncRNA_GIPSim,ddof=1)
LncRNA_SeqSim_mu=np.mean(LncRNA_SeqSim)#mean
LncRNA_SeqSim_sigma=np.std(LncRNA_SeqSim,ddof=1)

#plot
plt.figure(2)
plt.subplot(251)
disease_GIPSim_n, disease_GIPSim_bins, disease_GIPSim_patches = plt.hist(disease_GIPSim, num_bins, normed=True, color='steelblue')  #
disease_GIPSim_y = mlab.normpdf(disease_GIPSim_bins, disease_GIPSim_mu, disease_GIPSim_sigma)
plt.plot(disease_GIPSim_bins, disease_GIPSim_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of disease_GIPSim')
#设置两个坐标轴的范围  
plt.ylim(0,80)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(252)
diseasesim_DOSEwang_n, diseasesim_DOSEwang_bins, diseasesim_DOSEwang_patches = plt.hist(diseasesim_DOSEwang, num_bins, normed=True, color='blanchedalmond') 
diseasesim_DOSEwang_y = mlab.normpdf(diseasesim_DOSEwang_bins, diseasesim_DOSEwang_mu, diseasesim_DOSEwang_sigma)
plt.plot(diseasesim_DOSEwang_bins, diseasesim_DOSEwang_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of diseasesim_DOSEwang')
#设置两个坐标轴的范围  
plt.ylim(0,80)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(253)
diseasesim_funSim_n, diseasesim_funSim_bins, diseasesim_funSim_patches = plt.hist(diseasesim_funSim, num_bins, normed=True, color='lightgreen') 
diseasesim_funSim_y = mlab.normpdf(diseasesim_funSim_bins, diseasesim_funSim_mu, diseasesim_funSim_sigma)
plt.plot(diseasesim_funSim_bins, diseasesim_funSim_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of diseasesim_funSim')
#设置两个坐标轴的范围  
plt.ylim(0,80)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(254)
diseasesim_icodhprd_n, diseasesim_icodhprd_bins, diseasesim_icodhprd_patches = plt.hist(diseasesim_icodhprd, num_bins, normed=True, color='crimson') 
diseasesim_icodhprd_y = mlab.normpdf(diseasesim_icodhprd_bins, diseasesim_icodhprd_mu, diseasesim_icodhprd_sigma)
plt.plot(diseasesim_icodhprd_bins, diseasesim_icodhprd_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of diseasesim_icodhprd')
#设置两个坐标轴的范围  
plt.ylim(0,80)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(255)
diseasesim_jaccard_n, diseasesim_jaccard_bins, diseasesim_jaccard_patches = plt.hist(diseasesim_jaccard, num_bins, normed=True, color='deeppink') 
diseasesim_jaccard_y = mlab.normpdf(diseasesim_jaccard_bins, diseasesim_jaccard_mu, diseasesim_jaccard_sigma)
plt.plot(diseasesim_jaccard_bins, diseasesim_jaccard_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of diseasesim_jaccard')
#设置两个坐标轴的范围  
plt.ylim(0,80)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(256)
lncRNA_GIPSim_n, lncRNA_GIPSim_bins, lncRNA_GIPSim_patches = plt.hist(lncRNA_GIPSim, num_bins, normed=True, color='darkgrey') 
lncRNA_GIPSim_y = mlab.normpdf(lncRNA_GIPSim_bins, lncRNA_GIPSim_mu, lncRNA_GIPSim_sigma)
plt.plot(lncRNA_GIPSim_bins, lncRNA_GIPSim_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of lncRNA_GIPSim')
#设置两个坐标轴的范围  
plt.ylim(0,30)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

plt.subplot(258)
LncRNA_SeqSim_n, LncRNA_SeqSim_bins, LncRNA_SeqSim_patches = plt.hist(LncRNA_SeqSim, num_bins, normed=True, color='salmon') 
LncRNA_SeqSim_y = mlab.normpdf(LncRNA_SeqSim_bins, LncRNA_SeqSim_mu, LncRNA_SeqSim_sigma)
plt.plot(LncRNA_SeqSim_bins, LncRNA_SeqSim_y, 'r--')
plt.xlabel('Similarity')
plt.ylabel('probability')
plt.title(r'The distribution of LncRNA_SeqSim')
#设置两个坐标轴的范围  
plt.ylim(0,30)
plt.xlim(0,1)
# 设置图的底边距  
plt.subplots_adjust(bottom=0.15)

# the histogram of the data
#disease_GIPSim_n, disease_GIPSim_bins, disease_GIPSim_patches = plt.hist(disease_GIPSim, num_bins, normed=True, color='steelblue')  #
#diseasesim_DOSEwang_n, diseasesim_DOSEwang_bins, diseasesim_DOSEwang_patches = plt.hist(diseasesim_DOSEwang, num_bins, normed=True, color='blanchedalmond') 
#diseasesim_funSim_n, diseasesim_funSim_bins, diseasesim_funSim_patches = plt.hist(diseasesim_funSim, num_bins, normed=True, color='chocolate') 
#diseasesim_icodhprd_n, diseasesim_icodhprd_bins, diseasesim_icodhprd_patches = plt.hist(diseasesim_icodhprd, num_bins, normed=True, color='crimson') 
#diseasesim_jaccard_n, diseasesim_jaccard_bins, diseasesim_jaccard_patches = plt.hist(diseasesim_jaccard, num_bins, normed=True, color='firebrick') 
#lncRNA_GIPSim_n, lncRNA_GIPSim_bins, lncRNA_GIPSim_patches = plt.hist(lncRNA_GIPSim, num_bins, normed=True, color='lemonchiffon') 
#LncRNA_SeqSim_n, LncRNA_SeqSim_bins, LncRNA_SeqSim_patches = plt.hist(LncRNA_SeqSim, num_bins, normed=True, color='moccasin') 

plt.show()