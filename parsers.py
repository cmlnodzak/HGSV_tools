#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 10:15:25 2019

@author: cmlnodzak
"""


class Fasta:
    '''
    ====================================================================
                An example docstring for BINF Programming II
    Fasta:  a class for generating Fasta objects from microbiome dataset.
    ====================================================================
    Parameters: 
        Leave empty if one or more blank Fasta objects are required.
        Useful for parsing multifasta files.
        Alternatively, provide a single fasta as infile.
    ====================================================================
    ====================================================================
    Returns:
        A object of class Fasta which has an accession, rank, sequence, length
        and x - y coordinates. 
        
        Optionally, the user may run obj.setGC() on a Fasta object, obj, 
        to calculate a float of the percent GC content of the object.

    '''
    def __init__(self, infile=None):
        self.infile = infile
        self.accession = ""
        self.rank = ""
        self.sequence = ""
        self.length = 0
        self.x = 0.0
        self.y = 0.0
        self.GCcontent = None
        if self.infile != None:
            with open(self.infile, 'r') as fa:
                for line in fa.readlines():
                    if not line.startswith('>'):
                        self.sequence = self.sequence + line.rstrip()
                    elif line.startswith('>'):
                        self.accession = line.split()[0].lstrip('>')
                        self.rank = line.split()[1].split('=')[1]
                        self.x = line.split()[2].split('=')[1]
                        self.y = line.split()[3].split('=')[1]
                        self.length = line.split()[4].split('=')[1]    
    def getSeq(self):
        return self.sequence
    
    def setSeq(self,seq):
        self.sequence = self.sequence + line.rstrip()
    
    def getRank(self):
        return self.rank
    
    def setRank(self,line):
        self.rank = line.split()[1].split('=')[1]
    
    def getLength(self):
        return self.length
        
    def setLength(self,line):
        self.Length = line.split()[4].split('=')[1]    
    
    def getAccession(self):
        return self.accession
    
    def setAccession(self,line):
        self.accession = line.split()[0].lstrip('>')
    
    def getX(self):
        return self.x
    
    def setX(self,line):
        self.x = line.split()[2].split('=')[1]
    
    def getY(self):
        return self.y
    
    def setY(self,line):
        self.y = line.split()[3].split('=')[1]
    
    def setGC(self):
        if self.GCcontent == None:
            print("Running gcCalculator on " + self.accession)
            self.GCcontent = gcCalculator(self.sequence)
        return self.GCcontent
            
    def __repr__(self):
        return "A fasta object for: " + self.accession + " at X = " + self.x + " and Y = " + self.y
    
fastas = []

with open("BINF6112.lab.multi.fa", 'r') as mf:
    for line in mf.readlines():
        line = line.rstrip()
        if line.startswith(">"):
            faObj = Fasta()
            fastas.append(faObj)
            if faObj.sequence == "":
                faObj.setAccession(line)
                faObj.setRank(line)
                faObj.setLength(line)
                faObj.setX(line)
                faObj.setY(line)
                
        else:
            faObj.setSeq(line)
        
def gcCalculator(seq):
        G = seq.upper().count('G')
        C = seq.upper().count('C')
        percentGC = ((G + C) / len(seq))*100
        GCcontent = "%.2f" % round(percentGC,2)
        return GCcontent
       
        
for f in fastas:
    x = f.setGC()
    print(x)
    
    

        
         