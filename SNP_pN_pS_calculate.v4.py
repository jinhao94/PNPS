#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

import os,sys,argparse,shutil
def code_help():
    
    # args = sys.argv
    script_path = os.path.abspath(sys.argv[0])
    parser = argparse.ArgumentParser(description='A pipeline for calculating the PN/PS for genome.')
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--ffn", required = True, type=str, help="ffn path from Prodigal")
    parser.add_argument("-b", "--vcf", required = True, type=str, help="vcf file from bcftools")
    parser.add_argument("-o", "--outdir", required = True, type=str, help="name of outdir")
    parser.add_argument("-d", "--minDepth", type=int, help="minimum depth, default 20", default=20)
    parser.add_argument("-c", "--SNPcoverage", type=float,help="SNP coverage, default 0.05", default=0.05)
    # parser.add_argument("-ow", "--overwrite", default= False, action='store_true', help="Force overwrite the outdir")
    args = parser.parse_args()

    # if not args.ffn or not args.vcf or not args.outdir:
        # os.system('python %s -h'%script_path)
        # sys.exit()

    # return args.ffn,args.vcf,args.outdir,args.minDepth,args.SNPcoverage, args.overwrite
    return args.ffn,args.vcf,args.outdir,args.minDepth,args.SNPcoverage

ffnPath, vcfPath, outdir, minDepth, SNPcoverage  = code_help()

#================================================================================================

if not os.path.exists(outdir):
    os.makedirs(outdir)

import re
import pandas as pd
from math import ceil
import numpy as np

#ffnPath='../seq/4.Prodigal/C001'
#vcfPath='03.pN_pS/C001/filter.vcf'
#outdir='03.pN_pS/C001'
#minDepth=40
#SNPcoverage=0.05

#Data Preparation
#codon11Path='/share/data5/guorc/database/nt_20190330/codon.11_20191031'
codonNSPath='/nvmessdnode3/opt/.method/Genetic_Codes/sn_tot.list1'
codonNStranPath='/nvmessdnode3/opt/.method/Genetic_Codes/sn_tot.list2'

class count_pN_pS:
    
    def openf(self,Path,Line=False):
        with open(Path) as f:
            if Line:
                self.inputFile=f.readlines()
            else:
                self.inputFile=f.read()
    
    def ffnTreat(self,fPath):
        self.openf(fPath)
        self.ffn = self.inputFile.split('>')[1:]
        dic = {}
        for x in self.ffn:
            Tmp = x.split('\n',1)
            Seq = Tmp[1].replace('\n','')
            Tmp = Tmp[0].split(' ')
            ContigID = Tmp[0].rsplit('_',1)[0]
            if ContigID in dic:
                dic[ContigID].append({'SeqID':Tmp[0],'Start':int(Tmp[2]),\
                'End':int(Tmp[4]),'Direction':Tmp[6],'Seq':Seq})
            else:
                dic[ContigID] = [{'SeqID':Tmp[0],'Start':int(Tmp[2]),\
                'End':int(Tmp[4]),'Direction':Tmp[6],'Seq':Seq}]
        self.ffn = dic
    
    def vcfTreat(self,vPath,minD=False,SNPcov=False): ## filter vcf file by SNPcov and minD
        #self.vcf = pd.read_csv(vPath,sep='\t',header=None,comment='#')
        #self.vcf = self.vcf[self.vcf[0].isin(list(self.ffn.keys()))]
        self.vcf = vPath[vPath[0].isin(list(self.ffn.keys()))]
        self.I16 = [list(map(int,list(re.search('I16=(\d+),(\d+),(\d+),(\d+)',x).groups()))) for x in self.vcf[7]]
        Filter=[]
        for x in self.I16:
            if sum(x)>=minD and sum(x[2:])/sum(x)>=SNPcov:
                Filter.append(False)
            else:
                Filter.append(True)
        # print (self.vcf.iloc[Filter,:])
        # tmp_var = self.vcf.loc[Filter,'3']
        # self.vcf.loc[Filter,'4'] = tmp_var
        self.vcf.iloc[Filter,4]=self.vcf.iloc[Filter,3]
        # print (self.vcf.iloc[Filter,[3,4]])
    
    def Complement(self,Base):
        Base_C = {'A':'T','T':'A','G':'C','C':'G'}
        self.Base = Base_C[Base]
    
    def getSNPcodon(self):
        vcf_f=self.vcf[[0,1,3,4]].reset_index(drop=True)
        vcf_f.columns = ['ContigID','POS','Ref','Aln']
        Refcodon=[]
        SNPcodon=[]
        SeqID=[]
        Ncodon=[]
        Spos=[]
        Epos=[]
        self.SNPcodon = pd.DataFrame()
        for x,y in self.ffn.items():
            tab = vcf_f.loc[vcf_f['ContigID']==x]
            if tab.shape[0]==0:
                continue
            for z in y:
                tabf = tab[(tab['POS']>=z['Start']) & (tab['POS']<=z['End'])].reset_index(drop=True)
                if tabf.shape[0]==0:
                    continue
                
                Ref_codon = re.findall(r'.{3}',z['Seq'])
                SNP_seq = list(z['Seq'])
                POS = [i-z['Start']+1 for i in tabf['POS']]
                Ncodon = [ceil(j/3) for j in POS]
                # print ("aaa")
                if z['Direction']=='1':
                    Spos = [z['Start']+3*(j-1) for j in Ncodon]
                    Epos = [j+2 for j in Spos]
                    for i in range(len(POS)):
                        SNP_seq[POS[i]-1] = tabf.loc[i,'Aln'][0]
                    SNP_codon = re.findall(r'.{3}',''.join(SNP_seq))
                    SNP_codon = [SNP_codon[x-1] for x in Ncodon]
                    Ref_codon = [Ref_codon[x-1] for x in Ncodon] 
                else:
                    Epos = [z['Start']+3*(j-1) for j in Ncodon]
                    Spos = [j+2 for j in Epos]
                    for i in range(len(POS)):
                        self.Complement(tabf.loc[i,'Aln'][0])
                        SNP_seq[-POS[i]] = self.Base
                    SNP_codon = re.findall(r'.{3}',''.join(SNP_seq))
                    SNP_codon = [SNP_codon[-x] for x in Ncodon]
                    Ref_codon = [Ref_codon[-x] for x in Ncodon]
                    Ncodon = [int(len(z['Seq'])/3-(x-1)) for x in Ncodon]
                
                Tab = pd.DataFrame({'#SeqID':z['SeqID'],'Refcodon':Ref_codon,'SNPcodon':SNP_codon,'Ncodon':Ncodon,'Start':Spos,'End':Epos})
                Filter = Tab['Ncodon'].value_counts()
                Tab = Tab[Tab['Ncodon'].isin(Filter.index[Filter.values==3])].drop_duplicates()
                Tab['Ncodon'] = Tab['Ncodon']-1
                self.SNPcodon = pd.concat([self.SNPcodon, Tab], axis=0)
    
    def getpNpS(self,bPath,tPath):
        codon_NS = pd.read_csv(bPath,sep='\t',header=0)
        codon_NStran = pd.read_csv(tPath,sep='\t',header=0,usecols=[0,2,3,4])
        codon_NS.columns = ['Refcodon','N','S']
        codon_NStran.columns = ['Refcodon','SNPcodon','Nd','Sd']
        codon_NS['Refcodon'] = [x[0:3] for x in codon_NS['Refcodon']]
        codon_NStran['Refcodon'] = [x[0:3] for x in codon_NStran['Refcodon']]
        codon_NStran['SNPcodon'] = [x[0:3] for x in codon_NStran['SNPcodon']]
        self.SNPcodon = pd.merge(self.SNPcodon,codon_NS,how='left',on='Refcodon')
        self.SNPcodon = pd.merge(self.SNPcodon,codon_NStran,how='left',on=['Refcodon','SNPcodon'])
        self.SNPcodon = self.SNPcodon.fillna(0)
        if((sum(self.SNPcodon['Nd'])+sum(self.SNPcodon['Sd']))!=0):
            self.pN_pS = (sum(self.SNPcodon['Nd'])/sum(self.SNPcodon['N']))/(sum(self.SNPcodon['Sd'])/sum(self.SNPcodon['S']))
            print (sum(self.SNPcodon['Nd'],), sum(self.SNPcodon['N']), sum(self.SNPcodon['Sd']), sum(self.SNPcodon['S']))
        else:
            print (sum(self.SNPcodon['Nd']), sum(self.SNPcodon['N']), sum(self.SNPcodon['Sd']), sum(self.SNPcodon['S']))
            self.pN_pS = np.nan


##################################################################
absPath=os.path.abspath(ffnPath)
vPath=pd.read_csv(vcfPath,sep='\t',header=None,dtype=object)

Fun = count_pN_pS()
Group = vPath.groupby([0])
for group,dataframe in Group:
    ffn = absPath+'/'+group+'.ffn'
    vcf = dataframe.drop([0],axis=1).reset_index(drop=True)
    vcf.columns = range(vcf.shape[1])
    vcf[1] = vcf[1].astype("int")
    Fun.ffnTreat(ffn)
    Fun.vcfTreat(vcf,minDepth,SNPcoverage)
    Fun.getSNPcodon()
    Fun.getpNpS(codonNSPath,codonNStranPath)
    if Fun.SNPcodon.shape[0]==0:
        print('\n'+outdir+'/'+group+'.pnps'+': No Gene SNP!!!!\n')
    with open(outdir+'/'+group+'.pnps','w') as f:
        f.write('#pN/pS = '+str(Fun.pN_pS)+'\n')
    Fun.SNPcodon.to_csv(outdir+'/'+group+'.pnps',sep='\t',index=False, mode='a')


