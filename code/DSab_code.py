# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:18:37 2017

@author: bugatti
"""


import re
import swalign

#smith-waterman local alignment
def sw_one(query,refseq):
    match = 5
    mismatch = -4
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring,gap_penalty = -30,gap_extension_penalty = -1)
    alignment = sw.align(refseq, query)
    #score = alignment.score
    q_pos = alignment.q_pos
    q_end = alignment.q_end
    r_pos = alignment.r_pos
    #print q_pos, q_end, r_pos, r_end
    q_len = q_end-q_pos
    middle_q = q_pos+0.5*q_len
    middle_r = r_pos+0.5*q_len
    #print query,refseq
    #print middle_q, middle_r
    return middle_q, middle_r
    
#split query into motifs
def split_query(k, query):
    li_query = []
    li_q_loc = []
    for i in range(len(query)-k+1):
        li_query.append(query[i:i+k])
        loc = i+k*0.5+0.5
        li_q_loc.append(loc)
    return li_query, li_q_loc

#search and locate hotspots
def get_RGYW(ref, p, q):
    dic_hot = []
    iterator1 = re.finditer(p,ref)
    iterator2 = re.finditer(q,ref)
    for item in iterator1:
        loc = item.span()[0]+3
        dic_hot.append(loc)
    for item in iterator2:
        loc = item.span()[0]+2
        dic_hot.append(loc)
    dic_hot.sort()
    return dic_hot
    
#search and locate coldspots
def get_SYC(ref, p, q):
    dic_cold = []
    iterator1 = re.finditer(p,ref)
    iterator2 = re.finditer(q,ref)
    for item in iterator1:
        loc = item.span()[0]+3
        dic_cold.append(loc)
    for item in iterator2:
        loc = item.span()[0]+1
        dic_cold.append(loc)
    dic_cold.sort()
    return dic_cold

#read D region sequences
def read_db(path):
    dic_db = {}
    f = open(path, 'r')
    while 1:
        line = f.readline().strip()
        if line == '':
            break
        if line[0] == '>':
            fa = f.readline().strip()
            fa = fa.upper()
            #dic_db[li[1]] = fa
            dic_db[line[1:]] = fa
    f.close()
    return dic_db

#process db sequences
def db_process(dic_db, k, p, q, pp, qq):
    dic_pro_db = {}
    for item in dic_db:
        dic_pro_db[item] = []
        dic_pro_db[item].append(dic_db[item])
        dic_hot = get_RGYW(dic_db[item], p, q)
        dic_cold = get_SYC(dic_db[item], pp, qq)
        dic_pro_db[item].append(dic_hot)
        li_ref_motif, li_ref_loc = split_query(k, dic_db[item])
        dic_pro_db[item].append(dic_cold)
        dic_pro_db[item].append(li_ref_motif)
        dic_pro_db[item].append(li_ref_loc)
        #print dic_pro_db
    return dic_pro_db   #dic[name: seq, hotspot, coldspot, k-mers, location]

def motif_dis(mid_q, mid_r, loc_q, loc_r ):
    if loc_r <= mid_r:
        if loc_q <= mid_q:
            tmp_dis = abs((mid_r-loc_r)-(mid_q-loc_r))
        else:
            tmp_dis = (mid_r-loc_r)+(loc_r-mid_q)
    else:
        if loc_r <= mid_q:
            tmp_dis = (mid_q-loc_r)+(loc_r-mid_r)
        else:
            tmp_dis = abs((loc_r-mid_r)-(loc_r-mid_q))
    return tmp_dis
    
#map germline motif to the query
def best_match(li_pro_db, ger_num, mid_q, mid_r, li_query, li_q_loc, k, mismatch, tmp_max):
    li_score = []
    for i in range(len(li_query)):
        li_mut = []
        score = 0
        for j in range(len(li_query[i])):
            if li_query[i][j] == li_pro_db[3][ger_num][j]:
                pass
            else:
                loc_ger = int(li_pro_db[4][ger_num]-(k-1)/2)+j
                li_mut.append(loc_ger)
        if len(li_mut) <= mismatch:
            tmp_score = k-mismatch
            for m in range(len(li_mut)):
                if m in li_pro_db[1]:
                    tmp_score = tmp_score*0.8
                elif m in li_pro_db[2]:
                    tmp_score = tmp_score*0.2
                else:
                    tmp_score = tmp_score*0.5
            dis = motif_dis(mid_q, mid_r, li_q_loc[i], li_pro_db[4][ger_num])
            tmp_score = (1-dis*1.0/tmp_max)*tmp_score
            if tmp_score >= score:
                score = tmp_score
        li_score.append(score)
    best_score = max(li_score)
    #print best_score
    return best_score
        

#search the best reference sequence
def db_search(li_query, li_q_loc, dic_pro_db, k, mismatch, query, tmp_d):
    li_name = []
    li_score = []
    for item in dic_pro_db:
        item_score = 0
        ####similarity####
        mid_q, mid_r = sw_one(query, dic_pro_db[item][0])
        for i in range(len(dic_pro_db[item][3])):
            tmp_max = max(mid_q+len(dic_pro_db[item][0])-mid_r, mid_r+len(query)-mid_q)
            best_score = best_match(dic_pro_db[item], i, mid_q, mid_r, li_query, li_q_loc, k, mismatch, tmp_max)
            item_score += best_score
        li_name.append(item)
        li_score.append(item_score)
    #print li_name, li_score
    best_n = li_score.index(max(li_score))
    best_ref = li_name[best_n]
    if li_score[li_name.index(tmp_d)] == max(li_score):
        best_ref = tmp_d
    #print best_ref,max(li_score),tmp_d,li_score[li_name.index(tmp_d)]
    if max(li_score) == 0:
        best_ref = 'NA'
    #print best_ref, tmp_d
    return best_ref
            
            
################################################################################

def main_code(k, mismatch):
    print 'k',k,'mismatch',mismatch

    #hotspot pattern
    p = '[AT][GA]C[CT]'
    q = '[GA]G[CT][AT]'

    pp = '[CG][CT]C'
    qq = 'G[GA][CG]'

    #read D region sequences
    db_path = r'C:\D\VDJ\human_germline_for_igBlast\my_IGHD.fasta'
    dic_db = read_db(db_path)
    dic_pro_db = db_process(dic_db, k, p, q, pp, qq)

    #read query file
    f = open(r'C:\D\VDJ\99_fault_seq.fa','r')
    n = 0
    m = 0
    while 1:
        m += 1
        if m%500 == 0:
            print m, n
        line = f.readline().strip()
        if line == '':
            break
        li = line.split('\t')
        if len(li) < 4:
            continue
        query = li[3]
        tmp_d = li[0].split('|')[2]
        query = query.upper()
        li_query, li_q_loc = split_query(k, query)
        if len(li_query) == 0:
            continue
        best_ref = db_search(li_query, li_q_loc, dic_pro_db, k, mismatch, query, tmp_d)
        if best_ref == tmp_d:
            n += 1
    print n
    f.close()
    return n



#define the length and mismatch of motifs
g = open('C:\D\VDJ\99_fault_seq_D.txt','w')
for k in range(5, 8):
    for mismatch in range(0, 4):
        print k, mismatch
        g.write('k:'+str(k)+' mismatch:'+str(mismatch))
        g.write('\n')
        n = main_code(k, mismatch)
        g.write(str(n))
        g.write('\n')
        g.write('\n')
g.close()
print 'D!'