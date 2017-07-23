import numpy as np
from scipy import interp
from sklearn.metrics import roc_curve
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from sklearn.cross_validation import StratifiedKFold
import types
import scipy.io as sio

def plott(lab,score):
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, 100)
    for i, (train, test) in enumerate(cv):
    # Compute ROC curve and area the curve
    #通过roc_curve()函数，求出fpr和tpr，以及阈值
    #print (y[test])
    #print ( pred[test])
        fpr, tpr, thresholds = roc_curve(lab[test], score[test])
    
        #print (fpr)
        #print (tpr)
        #print (thresholds)
        mean_tpr += interp(mean_fpr, fpr, tpr)			#对mean_tpr在mean_fpr处进行插值，通过scipy包调用interp()函数
        mean_tpr[0] = 0.0 								#初始处为0
        roc_auc = auc(fpr, tpr)
        #画图，只需要plt.plot(fpr,tpr),变量roc_auc只是记录auc的值，通过auc()函数能计算出来
        #plt.plot(fpr, tpr, lw=1, label='ROC fold %d (area = %0.4f)' % (i, roc_auc))
    #画对角线
    plt.plot([0, 1], [0, 1], '--', color=(0.6, 0.6, 0.6))#, label='Luck'
    mean_tpr /= len(cv) 					#在mean_fpr100个点，每个点处插值插值多次取平均
    mean_tpr[-1] = 1.0 						#坐标最后一个点为（1,1）
    mean_auc = auc(mean_fpr, mean_tpr)		#计算平均AUC值
    #画平均ROC曲线
    #print mean_fpr,len(mean_fpr)
    #print mean_tpr
    plt.plot(mean_fpr, mean_tpr, 'k-',label='XXXX : AUC = %0.4f' % mean_auc, lw=1)
    #plt.xlim([-0.05, 1.05])
    #plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.show()

Result=sio.loadmat('result.mat')
print(Result)
DN_O_Sim=Result['nomarlized_DN_O_Sim_confidence_score_matrix']
DN_Sim=Result['nomarlized_DN_Sim_confidence_score_matrix']
lncN_O_Sim=Result['nomarlized_lncN_O_Sim_confidence_score_matrix']
lncN_Sim=Result['nomarlized_lncN_Sim_confidence_score_matrix']
lncRNA_Disease=Result['lncRNA_Disease_Matrix']

#DN_O=sio.loadmat('nomarlized_DN_O_Sim_confidence_score_matrix.mat')
#DN=sio.loadmat('nomarlized_DN_Sim_confidence_score_matrix.mat')
#lncN_O=sio.loadmat('nomarlized_lncN_O_Sim_confidence_score_matrix.mat')
#lncN=sio.loadmat('nomarlized_lncN_Sim_confidence_score_matrix.mat')
#lncRNA_Disease=sio.loadmat('lncRNA_Disease_Matrix.mat')

#DN_O_Sim=DN_O['DN_O_Sim']
#DN_Sim=DN_O['DN_Sim']
#lncN_O_Sim=lncN_O['lncN_O_Sim']
#lncN_Sim=lncN_O['lncN_Sim']
#lncRNA_Disease=lncRNA_Disease['lncRNA_Disease_Matrix']

#print ("DN_O_Sim")
#print (DN_O_Sim)
#print ("DN_Sim")
#print (DN_Sim)
#print ("lncN_O_Sim")
#print (lncN_O_Sim)
#print ("lncN_Sim")
#print (lncN_Sim)
#print ("lncRNA_Disease")
#print (lncRNA_Disease)

y=np.array(lncRNA_Disease)

DN_O_Sim_pred = np.array(DN_O_Sim)
DN_Sim_pred = np.array(DN_Sim)
lncN_O_Sim_pred = np.array(lncN_O_Sim)
lncN_Sim_pred = np.array(lncN_Sim)

y=y.reshape((1,18603));
y=y[0];
#print (y)

DN_O_Sim_pred=DN_O_Sim_pred.reshape((1,18603));
DN_O_Sim_pred=DN_O_Sim_pred[0];
DN_Sim_pred=DN_Sim_pred.reshape((1,18603));
DN_Sim_pred=DN_Sim_pred[0];
lncN_O_Sim_pred=lncN_O_Sim_pred.reshape((1,18603));
lncN_O_Sim_pred=lncN_O_Sim_pred[0];
lncN_Sim_pred=lncN_Sim_pred.reshape((1,18603));
lncN_Sim_pred=lncN_Sim_pred[0];

#print (pred)

cv = StratifiedKFold(y, n_folds=10)

plott(y,DN_O_Sim_pred)
#plott(y,DN_Sim_pred)
#plott(y,lncN_O_Sim_pred)
#plott(y,lncN_Sim_pred)
