class Config(object):
    def __init__(self):
        ## graph data
        ##self.file_path = "GraphData/blogCatalog3.txt"
        ##self.file_path = "GraphData/DlncN117.txt"
        self.file_path = "GraphData/lncDN159.txt"
        ##self.label_file_path = "GraphData/blogCatalog-groups.txt"
        ##self.label_file_path = "GraphData/DlncN117-groups.txt"
        #self.label_file_path = "GraphData/lncDN159-groups.txt"
        self.label_file_path = "GraphData/lncDN159-groups-nonlabel.txt"
        
        ## embedding data
        ##self.embedding_filename = "embeddingResult/blogcatalog31" 
        ##self.embedding_filename = "embeddingResult/DlncN117"
        self.embedding_filename = "embeddingResult/lncDN159-embedding-nonlabel"
        ## hyperparameter
        self.struct = [None, 1000, 100]
        ## the loss func is  // gamma * L1 + alpha * L2 + reg * regularTerm // 
        self.alpha = 5
        self.gamma = 1
        self.reg = 1
        ## the weight balanced value to reconstruct non-zero element more.
        self.beta = 15
        
        ## para for training
        #self.batch_size = 1000
        #self.epochs_limit = 1000
        #self.learning_rate = 0.001
        self.batch_size = 1000
        self.epochs_limit = 3000
        self.learning_rate = 0.01
        
        self.DBN_init = False
        self.sparse_dot = False
        self.ng_sample_ratio = 0.1 # negative sample ratio
        
        
