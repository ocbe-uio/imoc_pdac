import itertools
from functools import reduce
import numpy as np
import pandas as pd
from gensim.models import Word2Vec


class RandomWalker:
    def __init__(self, graphs, nodes, verbose, random_state):
        """
        :param graphs: list of graph
        :param nodes: all nodes in the graphs
        :param verbose: whether print detail information
        """
        self.Graphs = graphs
        self.alias_nodes=[]
        self.nodes=nodes
        self.verbose=verbose
        if random_state is None:
            random_state = int(np.random.default_rng().integers(10000))
        self.random_state=random_state

    def walk(self, walk_length, start_node):
        Graphs = self.Graphs
        alias_nodes = self.alias_nodes
        walk = [start_node]
        while len(walk) < walk_length:
            cur = walk[-1]
            cur_nbrs = [list(G.neighbors(cur)) if(G.has_node(cur)) else [] for G in Graphs ]
            cand=[i for i, e in enumerate(cur_nbrs) if len(e) != 0]
            if(len(cand)>0):
                self.random_state = self.random_state + np.random.default_rng(self.random_state).integers(10000)
                select=np.random.default_rng(self.random_state).choice(cand)
                walk.append(
                    cur_nbrs[select][self._alias_sample(*alias_nodes[select][cur])])
            else:
                break
        return walk

    def simulate_walks(self ,num_walks ,walk_length ,workers=1):

        nodes=self.nodes

        results = [
            self._simulate_walks(nodes, num, walk_length) for num in
            RandomWalker._partition_num(num_walks, workers)]

        walks = list(itertools.chain(*results))

        return walks

    def _simulate_walks(self, nodes, num_walks, walk_length):
        walks = []
        for _ in range(num_walks):
            self.random_state = self.random_state + np.random.default_rng(self.random_state).integers(10000)
            nodes = np.random.default_rng(self.random_state).permutation(nodes)
            for v in nodes:
                walks.append(self.walk(
                    walk_length=walk_length, start_node=v))
        return walks

    def preprocess_transition_probs(self,G):
        """
        Preprocessing of transition probabilities for guiding the random walks.
        """
        alias_nodes = {}
        L=len(G.nodes())/100
        for i,node in enumerate(G.nodes()):
            if(self.verbose!=0 and i%L==0):
                print(i/len(G.nodes()))
            unnormalized_probs = [G[node][nbr].get('weight', 1.0)
                                  for nbr in G.neighbors(node)]
            norm_const = sum(unnormalized_probs)
            normalized_probs = [
                float(u_prob)/norm_const for u_prob in unnormalized_probs]
            alias_nodes[node] = RandomWalker._create_alias_table(normalized_probs)

        self.alias_nodes.append(alias_nodes)


    @staticmethod
    def _create_alias_table(area_ratio):
        """
        :param area_ratio: sum(area_ratio)=1
        :return: accept,alias
        """
        l = len(area_ratio)
        accept, alias = [0] * l, [0] * l
        small, large = [], []
        area_ratio_ = np.array(area_ratio) * l
        for i, prob in enumerate(area_ratio_):
            if prob < 1.0:
                small.append(i)
            else:
                large.append(i)

        while small and large:
            small_idx, large_idx = small.pop(), large.pop()
            accept[small_idx] = area_ratio_[small_idx]
            alias[small_idx] = large_idx
            area_ratio_[large_idx] = area_ratio_[large_idx] - \
                (1 - area_ratio_[small_idx])
            if area_ratio_[large_idx] < 1.0:
                small.append(large_idx)
            else:
                large.append(large_idx)

        while large:
            large_idx = large.pop()
            accept[large_idx] = 1
        while small:
            small_idx = small.pop()
            accept[small_idx] = 1

        return accept, alias


    def _alias_sample(self, accept, alias):
        """
        :param accept:
        :param alias:
        :return: sample index
        """
        N = len(accept)
        self.random_state = self.random_state + np.random.default_rng(self.random_state).integers(10000)
        i = int(np.random.default_rng(self.random_state).random()*N)
        self.random_state = self.random_state + np.random.default_rng(self.random_state).integers(10000)
        r = np.random.default_rng(self.random_state).random()
        if r < accept[i]:
            return i
        else:
            return alias[i]


    @staticmethod
    def _partition_num(num, workers):
        if num % workers == 0:
            return [num // workers] * workers
        else:
            return [num // workers] * workers + [num % workers]


class Embedding:

    def __init__(self, graphs, workers=1,verbose=0, random_state= None):

        self.random_state = random_state
        self.graphs = graphs
        self._embeddings = {}
        self.nodes = list(reduce(lambda x, y: x | y, [set(g.nodes()) for g in self.graphs]))
        self.walker = RandomWalker(graphs,nodes=self.nodes,verbose=verbose, random_state=random_state)
        self.workers=workers
        self.verbose=verbose
        for g in graphs:
            self.walker.preprocess_transition_probs(g)

    def train(self, embed_size=128, window_size=5,epoch=5, **kwargs):

        kwargs["sentences"] = self.sentences
        kwargs["min_count"] = kwargs.get("min_count", 0)
        kwargs["vector_size"] = embed_size
        kwargs["sg"] = 1  #select sg mode
        kwargs["hs"] = 0  #not use Hierarchical Softmax
        kwargs["workers"] = self.workers
        kwargs["window"] = window_size
        kwargs["epochs"] = epoch
        kwargs["seed"] = self.random_state

        model = Word2Vec(**kwargs)

        self.w2v_model = model

        return model

    def sample_sentence(self,num_walks,walk_length):
        #print("sampling sentences...")
        self.sentences = self.walker.simulate_walks(
            num_walks=num_walks, walk_length=walk_length, workers=self.workers)

    def get_embeddings(self,samples=None):
        if self.w2v_model is None:
            print("model not train")
            return np.array([])

        self._embeddings = {}
        if(samples is None):
            samples=self.nodes
        for word in samples:
            self._embeddings[word] = self.w2v_model.wv[word]

        embeddings=self._embeddings
        emb_list = []
        for k in embeddings:
            emb_list.append(embeddings[k])
        emb_list = np.array(emb_list)
        emb_list = emb_list / np.sqrt((emb_list * emb_list).sum(axis=1).reshape(-1, 1))
        return pd.DataFrame(emb_list,index=list(samples))
