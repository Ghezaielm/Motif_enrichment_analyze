import re
import numpy as np 
import matplotlib.pyplot as plt 
import os

class MEA(): 
    def __init__(self): 
        self.length = 10
        self.sequence = 0
        self.results = {}
        self.query_pfm = 0
        
    def getPfms(self): 
        self.motif_names =  os.listdir("Motifs")
        self.dics = {}
        for motif in self.motif_names:
            self.dic = {"A":[],"C":[],'G':[],"T":[]}
            with open(os.path.join("Motifs",motif)) as file:
                lines = file.readlines()
                for idx, nt in enumerate(self.dic):
                    line = lines[idx]
                    vals = [float(i) for i in re.findall("[0-9].[0-9]{1,8}",lines[idx])]
                    self.dic[nt]+=[i/max(vals) if max(vals)>0 else 0 for i in vals]
            self.dics[motif.split(".")[0]]=self.dic
            
    def getScore(self): 
        perfect_target = [i.index(max(i)) for i in self.pfm]
        perfect_target = [[1 if j==perfect_target[idx] else 0 for j in range(len(i))] for idx,i in enumerate(self.pfm)]
        worst_target = [i.index(min(i)) for i in self.pfm]
        worst_target = [[1 if j==worst_target[idx] else 0 for j in range(len(i))] for idx,i in enumerate(self.pfm)]
        trace_perfect = np.trace(np.dot(np.array(self.pfm),np.array(perfect_target).T))
        trace_worst = np.trace(np.dot(np.array(self.pfm),np.array(worst_target).T))
        trace_query  = np.trace(np.dot(np.array(self.pfm),np.array(self.query_pfm).T))
        sim = 1-((trace_perfect-trace_query-trace_worst)/(trace_perfect))
        return sim
        
    def getQueryPfm(self,target): 
        for id,motif in enumerate(self.dics):
            print(motif,id)
            curr_motif = self.dics[motif]
            motif_length = len(curr_motif["A"])
            self.pfm = np.array([curr_motif[i] for i in curr_motif]).T.tolist()
            scores = []
            for win in range(0,len(target[0])-motif_length):
                self.sequence = target[0][win:win+motif_length]
                self.query_pfm = [[len(re.findall("A",target[0][idx])),
                        len(re.findall("C",self.sequence[idx])),
                        len(re.findall("G",self.sequence[idx])),
                        len(re.findall("T",self.sequence[idx]))]
                        for idx in range(len(self.sequence))]

                scores.append(MEA.getScore(self))
            self.results[motif]=np.mean(scores)/np.std(scores)
    
              
o = MEA()
o.getPfms()
target = ["GAGCCCCGAGCCCCGAGCCCC"]
o.getQueryPfm(target)
#o.generatePfm()
