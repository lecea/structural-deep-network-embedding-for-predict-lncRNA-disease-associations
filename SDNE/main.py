'''
Reference implementation of SDNE

Author: Xuanrong Yao, Daixin wang

for more detail, refer to the paper:
SDNE : structral deep network embedding
Wang, Daixin and Cui, Peng and Zhu, Wenwu
Knowledge Discovery and Data Mining (KDD), 2016
'''

#!/usr/bin/python2
# -*- coding: utf-8 -*-



from config import Config
from graph import Graph
from model.sdne import SDNE
from utils.utils import *
import scipy.io as sio
import time
import scipy.io as sio

if __name__ == "__main__":
    config = Config()
    
    graph_data = Graph(config.file_path, config.ng_sample_ratio)
    graph_data.load_label_data(config.label_file_path)
    config.struct[0] = graph_data.N
    
    model = SDNE(config)    
    model.do_variables_init(config.DBN_init)

    last_loss = np.inf
    converge_count = 0
    time_consumed = 0
    epochs = 0
    batch_n = 0
    while (True):
        mini_batch = graph_data.sample(config.batch_size)
        st_time = time.time()
        model.fit(mini_batch)
        batch_n += 1
        time_consumed += time.time() - st_time
        #print ("Mini-batch : %d fit time: %.2f" % (batch_n, time_consumed))
        if graph_data.is_epoch_end:
            epochs += 1
            loss = 0
            embedding = None
            while (True):
                mini_batch = graph_data.sample(config.batch_size, do_shuffle = False)
                loss += model.get_loss(mini_batch)
                if embedding is None:
                    embedding = model.get_embedding(mini_batch)
                else:
                    embedding = np.vstack((embedding, model.get_embedding(mini_batch)))
                
                if graph_data.is_epoch_end:
                    break
            
            #print ("Epoch : %d Loss : %.3f, Train time_consumed : %.3fs" % (epochs, loss, time_consumed))
            if epochs % 100 == 0:
                ##check_link_reconstruction(embedding, graph_data, [10000,30000,50000,70000,90000,100000])
                check_link_reconstruction(embedding, graph_data, [50,100,150,200,250,300])#[200,400,600,800,1000,1200])
                data = graph_data.sample(graph_data.N, with_label = True)
                check_classification(model.get_embedding(data), data.label, test_ratio = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5])#[0.9, 0.7, 0.5, 0.3, 0.1])
            if (loss > last_loss):
                converge_count += 1
                if converge_count > 500:
                    print ("model converge terminating")
                    break
            if epochs > config.epochs_limit:
                print ("exceed epochs limit terminating")
                break
            last_loss = loss
        
    
    sio.savemat(config.embedding_filename + '_embedding.mat',{'embedding':embedding})
