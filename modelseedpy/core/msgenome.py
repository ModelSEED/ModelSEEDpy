import logging

import re
import copy
from cobra.core.dictlist import DictList

logger = logging.getLogger(__name__)

def load_fasta:
    


class GenomeFactory():
    def __init__(self,)
        self.features = DictList;
    
    
    #Creates a genome from a FASTA file - could be a protein feature fasta, a DNA feature fasta, or a contig fasta file 
    def build_from_fasta(self,filename,contigs = 0):
        sequences = load_fasta(filename)
        
        
        
        return genome
        
    def build_from_dna_fasta(self,filename):
            
    def build_from_kbase_gto:
        
    def build_from_kbase_id:
        
    def build_from_protein_sequence_hash:
    
    def build_from_protein_sequence_hash:
    
class MSFeature():
    def __init__(self,id,data):
        self.id = id
        self.length = None
        self.protein_sequence = None
        self.dna_sequence = None
        self.location = None
        if "location" in data:
            self.location = data["location"]
        if "dna_sequence" in data:
            self.dna_sequence = data["dna_sequence"]
            self.length = len(self.dna_sequence)/3
        if "protein_sequence" in data:
            self.protein_sequence = data["protein_sequence"]
            self.length = len(self.protein_sequence)
        elif self.dna_sequence != None:
            self.protein_sequence = self.translate()
            self.length = len(self.protein_sequence)
        

        
        
        for field in fields:
            if field in data:
                self.
    