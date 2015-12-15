import pysam as ps
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet
from matplotlib import pyplot as plt
from scipy import stats, optimize


def geneFrame(genbank_file):
    file = SeqIO.read(genbank_file,'gb')
    genes = []
    feature_type = []
    name = []
    product = []
    func = []
    strand = []
    start = []
    stop = []
    aaseq = []
    seq = []
    genome_seq_df = pd.DataFrame({'sequence':list(str(file.seq))},index=range(1,len(str(file.seq))+1))
    for feature in file.features:
        if (feature.type == 'CDS' or feature.type == 'rRNA' or feature.type == 'tRNA' or feature.type == 'ncRNA'):
            feature_type.append(feature.type)
            genes.append(feature.qualifiers['locus_tag'][0])
            name.append(feature.qualifiers['gene'][0])
            if 'product' in feature.qualifiers:
                product.append(feature.qualifiers['product'][0])
            else:
                product.append('N/A')
            if 'function' in feature.qualifiers:
                func.append(feature.qualifiers['function'][0])
            else:
                func.append("N/A")
            if 'translation' in feature.qualifiers:
                aaseq.append(feature.qualifiers['translation'][0])
            else:
                aaseq.append("N/A")
            if feature.strand == 1:
                strand.append("plus")
                start.append(feature.location.start.real+1)
#This means that we will have to use .ix when calling from pandas
                stop.append(feature.location.end.real)
                seq.append(str(file.seq[feature.location.start.real: feature.location.end.real]))
            elif feature.strand == -1:
                strand.append("minus")
                start.append(feature.location.start.real+1)
                stop.append(feature.location.end.real)
                seq.append(str(file.seq[feature.location.start.real: feature.location.end.real].reverse_complement()))
    df = pd.DataFrame({"gene": genes, "name": name, "type": feature_type, "product": product, "function": func, "strand": strand, "start": start, "stop": stop, "seq":seq, "aaseq": aaseq}, 
                                    columns = ["gene", "name", "type", "function", "product", "strand", "start", "stop", "seq", "aaseq"])
    df = df.set_index("gene")
    return df, genome_seq_df






import sys, time
try:
    from IPython.display import clear_output
    have_ipython = True
except ImportError:
    have_ipython = False

class ProgressBar:
    def __init__(self, iterations, secondary_iter=0):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 40
        self.secondary_iter = secondary_iter
        self.__update_amount(0)
        if have_ipython:
            self.animate = self.animate_ipython
        else:
            self.animate = self.animate_noipython

    def animate_ipython(self, iter):
        print '\r', self,
        sys.stdout.flush()
        self.update_iteration(iter + 1)

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        if self.secondary_iter>0:
            self.prog_bar += '  %d of %s complete (%d)' % (elapsed_iter, self.iterations, self.secondary_iter)
        else:
            self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)
