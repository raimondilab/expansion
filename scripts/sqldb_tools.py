
import os,sys,operator, math, argparse, pickle
sys.path.insert(0, "/projects/bioinformatics/users_scripts/")
from io import StringIO
import gzip, sqlite3, glob
import numpy as np
import glob, itertools
import timeit
from metric_tools import Overlap

fs="\t"

####A class to retrieve and manipulate interpro annotations from an in-house SQLlite database
class interprodb(object):

  def __init__(self):
    self.output={}

  ######Function to extract Intepro domain instances for a given input of Uniprot Acession 
  def extract(self, idslist):
    
    db="/projects/bioinformatics/DB/Interpro/interpro.db"
    conn = sqlite3.connect(db)

    c = conn.cursor()
    sql_query=""
    outstr=""
    for n in idslist:
      sql_query = 'select * from interpro where ac = '+"'"+ str(n)+"'"
      c.execute(sql_query)
      data=c.fetchall()
      if len(data) == 0:
        N=n.upper()
        sql_query = 'select * from interpro where ac = '+"'"+ str(N)+"'"
    
      for line in c.execute(sql_query):
        uniac=line[0]
        self.output[line[0]]=line[1]

    conn.close()
  
    return self.output

  def make_unique(self, overlap_cutoff=0.5):
    ip2type={}
    with open('/projects/bioinformatics/DB/Interpro/ip2type.pickle', 'rb') as handle:
      ip2type = pickle.load(handle)
    offset=2 ####Offset required to merge contiguous domains with the same Interpro IDs
    final_dom={}
    final_hier={}
    for uac in self.output.keys():
      domains={}
      dom_track={}
      for dom in self.output[uac].split("\n"):
        if len(dom.split(fs)) > 1:
          did=dom.split(fs)[0]
          dname=dom.split(fs)[1]
          d2ndid=dom.split(fs)[2]
          dstart=int(dom.split(fs)[3])
          dstop=int(dom.split(fs)[4].strip("\n"))
          if ip2type[did] == "Domain" or ip2type[did] == "Homologous_superfamily":
            if did not in dom_track:
              dom_track[did]=1
              did_lab=did+"_"+str(dom_track[did])
            else:
              dom_track[did]+=1
              did_lab=did+"_"+str(dom_track[did])
            #domains[did_lab]=[dname,dstart,dstop]
            ####If Interpro ID instance is not yet present, create it
            if dom_track[did] == 1:
              domains[did]=[dname,dstart,dstop]
            ###If Interpro ID instance is already present, check for the boundaries and eventually extend it
            ###Using an offset on the start and stop of each domain to incorporate flanking domains with the same Interpro ID
            elif dom_track[did] > 1:
              if dstart in range(domains[did][1]-offset, domains[did][2]+offset) or dstop in range(domains[did][1]-offset, domains[did][2]+offset):
                if dstart < domains[did][1]:
                  domains[did][1]=dstart
                if dstop > domains[did][2]:
                  domains[did][2]=dstop
              else:
                domains[did_lab]=[dname,dstart,dstop] 

      ###Check those redundant domain instances overlapping for more than an overlap cutoff (by default 50% with respect to the shortest range)
      ###The shorter domains instances overlapping for more than 50% with larger domains are discarded
      dom_list=list(domains.keys())
      blacklist=[]
      dom_hierarchy={}
      for ii in range(len(dom_list)):
        d1=dom_list[ii]
        d1_pre=d1.split("_")[0]
        n1=domains[d1][0]
        #if n1 not in dom_hierarchy:
        #  dom_hierarchy[n1]=n1
        for jj in range(ii+1, len(dom_list)):
           d2=dom_list[jj]
           d2_pre=d2.split("_")[0]
           n2=domains[d2][0]
           #if n2 not in dom_hierarchy:
           #  dom_hierarchy[n2]=n2
           range1=range(domains[d1][1], domains[d1][2])
           range2=range(domains[d2][1], domains[d2][2])
           O , order = Overlap(range1, range2)
           if O > overlap_cutoff:
             #print (d1, domains[d1], d2, domains[d2], O)
             if order == 0:
               blacklist.append(d1)
               if n2 not in dom_hierarchy:
                 dom_hierarchy[n2]=[n1]
               elif n1 not in dom_hierarchy[n2]:
                 dom_hierarchy[n2].append(n1)
             else:
               blacklist.append(d2)
               if n1 not in dom_hierarchy:
                 dom_hierarchy[n1]=[n2]
               elif n2 not in dom_hierarchy[n1]:
                 dom_hierarchy[n1].append(n2)
         ###This piece of code checks if two domains, with the same Interpro ID are contigous and then merges them by creating 
         
      final_dom[uac]={}
      final_hier[uac]={}
      
      for d in dom_list:
        if d not in blacklist:
          final_dom[uac][d]=domains[d]
          n=domains[d][0]
          if n in dom_hierarchy:
            for nn in dom_hierarchy[n]:
              final_hier[uac][nn]=n
          else:
            final_hier[uac][n]=n

      #for d in dom_hierarchy.keys():
      #  print (d, dom_hierarchy[d])
    return final_dom, final_hier




class pdb2sprotdb(object):

  def __init__(self):
    self.output={}

  ######Function to extract Intepro domain instances for a given input of Uniprot Acession 
  def extract(self, idslist):
    
    db="/projects/bioinformatics/DB/blastdb/pdb_sprot_mappings.db"
    conn = sqlite3.connect(db)

    c = conn.cursor()
    sql_query=""
    outstr=""
    for n in idslist:
      sql_query = 'select * from pdbsprot where uniac = '+"'"+ str(n)+"'"
      c.execute(sql_query)
      data=c.fetchall()
      if len(data) == 0:
        N=n.upper()
        sql_query = 'select * from pdbsprot where uniac = '+"'"+ str(N)+"'"

      outstr=""    
      for line in c.execute(sql_query):
        #self.output[line[0]]=[line[0], line[1], line[2]]
        #print (line)
        #outstr=+line[1]+"\t"+line[2]+"\t"+line[3]+"\n"
        pdb=line[0]
        chain=line[1].split("|")[0]
        res=line[1].split("|")[1]
        res_pos=line[1].split("|")[1][1:]
        res_aa=line[1].split("|")[1][1]
        uac=line[2]
        uac_pos=line[3]
          ######################pdbid########chain+pos#####uniac#######uniac_pos
          #self.output[line[0]]=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n"
        #elif line[0] in self.output:
        #  self.output[line[0]]+=line[0]+"\t"+line[1]+"\t"+line[2]+"\t"+line[3]+"\n"
        if pdb not in self.output:
          self.output[pdb]={}
          self.output[pdb][chain]={}
          self.output[pdb][chain][res_pos]=[res, uac, uac_pos]
        elif pdb in self.output and chain not in self.output[pdb]:
          self.output[pdb][chain]={}
          self.output[pdb][chain][res_pos]=[res, uac, uac_pos]
        elif pdb in self.output and chain in self.output[pdb]:
          self.output[pdb][chain][res_pos]=[res, uac, uac_pos]
        
          
    conn.close()
  
    return self.output






















#########################################################################################################################################################################################################################################################################################
##########This is the old wrapper function for legacy applications
def ExtractInterpro(db, idslist):
  output={}

  #fastaoutfile=open(fastaout,"w")

  conn = sqlite3.connect(db)

  c = conn.cursor()
  ###We need to retrieve all the pdb2uniprot matches for the input XL
  sql_query=""
  outstr=""
  #sql_query = 'select * from swissprot where ' + ' or '.join(("id_ac = '" + str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'" for n in idslist)) 
  for n in idslist:
    sql_query = 'select * from interpro where ac = '+"'"+ str(n)+"'"
    c.execute(sql_query)
    data=c.fetchall()
    if len(data) == 0:
      N=n.upper()
      sql_query = 'select * from interpro where ac = '+"'"+ str(N)+"'"
    
    for line in c.execute(sql_query):
      uniac=line[0]
      output[line[0]]=line[1]

  conn.close()
  
  return output

