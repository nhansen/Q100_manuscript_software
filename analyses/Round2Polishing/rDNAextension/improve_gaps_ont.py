#PUBLIC DOMAIN NOTICE
#
#This software is "United States Government Work" under the terms of the United
#States Copyright Act. It was written as part of the authors' official duties
#for the United States Government and thus cannot be copyrighted. This software
#is freely available to the public for use without a copyright
#notice. Restrictions cannot be placed on its present or future use.
#
#Although all reasonable efforts have been taken to ensure the accuracy and
#reliability of the software and associated data, the National Human Genome
#Research Institute (NHGRI), National Institutes of Health (NIH) and the
#U.S. Government do not and cannot warrant the performance or results that may
#be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#disclaim all warranties as to performance, merchantability or fitness for any
#particular purpose.
#
#Please cite Dmitry Antipov in any work or product based on this material.

#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
import graph_functions

#improve (rDNA) gaps with ONT alignments.
#required additional graphaligner run (otherwise only 1-smth alignments)

rukki_file = sys.argv[1]
ont_aligns = sys.argv[2]
graph = sys.argv[3]
#edges shorter are not distinctive
cutoff = int(sys.argv[4])
#edges with this and higher coverage are considered repeats
#can be automatized but do not want
avgcov = float(sys.argv[5])
#S       utig4-0 *       LN:i:131601635  RC:i:8039543882 ll:f:61.1
manual = sys.argv[6]
manual_repeat = set()
manual_unique = set()

#chromosomes to ignore
forbidden = set()
for line in open (manual):
    arr = line.strip().split()
    if arr[0] == "unique":
        for i in range (1, len(arr)):
            manual_unique.add(arr[i])
    elif arr[0] == "repeat":
        for i in range (1, len(arr)):
            manual_repeat.add(arr[i])
    elif arr[0] == "forbidden":
        for i in range (1, len(arr)):
            forbidden.add(arr[i])
    else:
        print("corrupted manual file")
        exit()
    
nodelens = {}
nodecov = {}
link_len = {}
for line in open (graph):
    arr = line.split()
    if arr[0] == "S":
        node = arr[1]
        lens = int(arr[3].split(':')[2])
        for pref in {'>', '<'}:
            nodelens[pref + node] = lens
            nodecov[pref + node] = float(arr[5].split(':')[2])
#L       utig4-984       +       utig4-1202      +       2929M           
    if arr[0] == "L":
        if arr[2] == "+":
            f1 = ">" + arr[1]
            b2 = "<" + arr[1]
        else:
            f1 = "<" + arr[1]
            b2 = ">" + arr[1]

        if arr[4] == "+":
            b1 = ">" + arr[3]
            f2 = "<" + arr[3]
        else:
            b1 = "<" + arr[3]
            f2 = ">" + arr[3]
        if not f1 in link_len:
            link_len[f1] = {}
        if not f2 in link_len:
            link_len[f2] = {}
    
        llen = 0
        if arr[4][:-1].isdigit():
            llen = int(arr[4][:-1])

        link_len[f1][b1] = llen
        link_len[f2][b2] = llen



def parse_path(rukki_path):
    global weights_map
    res = 0
    separators = ">|<|\[|\]|,"
    pattern = r'(<utig4-\d+)|(>utig4-\d+)|([.*])'
    edges = re.split(pattern, rukki_path)
    edges = list(filter(None, edges))
    return edges

reads = []
for line in open (ont_aligns):
    seq = line.strip().split()[5]
    
    arr = parse_path(seq)
    if len(arr) > 1:
        reads.append(seq)    
#        reads.append(arr.copy())    

        arr.reverse()
        for i in range (0, len(arr)):
            if arr[i][0] == '>':            
                arr[i] = '<' + arr[i][1:]
            elif arr[i][0] == '<':
                arr[i] = '>' + arr[i][1:]
    
        reads.append("".join(arr))
#    reads.append(arr.copy())    
sys.stderr.write("Parsed input\n")
directions = [-1, 1] #before and after gap
for line in open (rukki_file):
    arr = line.split()
    if len(arr) < 2:
        continue
    if arr[0] in forbidden:
        continue
    nodes = parse_path(arr[1])
    for i in range (0, len(nodes)):
        if nodes[i][0] == '[':
            print (f"found gap {nodes[i]}")
            sys.stderr.write(f"Closing gap for {arr[0]}\n")
#1 -extending from left to right, inside gap, -1 otherwise            
            addon = {}
            for direction in directions:
                print ()
                active_reads = reads.copy()
                pos = i + direction
                to_search = ""
                active_node_list = []                
                print (f"{i}  {pos}   {direction}  {len(nodes)}   ")
                while pos >= 0 and pos < len(nodes):
                    if nodes[pos] in nodelens:
                        if direction > 0:
                            to_search = to_search + nodes[pos]
                            active_node_list.append(nodes[pos])
                        else:
                            to_search = nodes[pos] + to_search
                            active_node_list = [nodes[pos]] + active_node_list
                        #cov small - likely repeat coverage drop. To be fixed with utig1-
                        if (nodes[pos][1:] in manual_unique) or (nodecov[nodes[pos]] > avgcov * 0.5 and nodelens[nodes[pos]] > cutoff and nodecov[nodes[pos]] < avgcov * 1.5 and (not (nodes[pos][1:] in manual_repeat))):
                            break
                    pos += direction  
                reuse = len(active_node_list)
                progress = True
                shift_to_uniq = 0

                
                while True:
                    if not (progress):
                        continue
                    progress = False
                    print (f"searching... {to_search}")
                    active = ""
                    for nd in active_node_list:
                        active += f"{nd}:{nodelens[nd]}:{nodecov[nd]}  " 
                    print (active)
                    useful_reads = []
                    str_reads = [] 
                    for read in reads:
                        if read.find(to_search)!= -1:
                            align = read
                            uuu = parse_path(align)
                            str_reads.append(read)
                            if len(uuu) > len(active_node_list): #otherwise no sence
                                useful_reads.append(uuu) 
                    
                    extensions = {}
                    stopped = 0
#                    for useful_read in useful_reads:#useful_reads:                    
                    for useful_read in useful_reads:#useful_reads:                    
 #                       useful_read = parse_path(useful_read_str)
                        for ii in range (0, len(useful_read)):
                            found = True
                            for k in range (shift_to_uniq, len(active_node_list)):
                                if (ii + k >= len(useful_read) + shift_to_uniq or active_node_list[k] != useful_read[ii + k - shift_to_uniq]):
                                    found = False
                                    break
                            if found:
                                if direction == -1:
                                    next_pos = ii + len(active_node_list) - shift_to_uniq
                                else:
                                    next_pos = ii - 1 - shift_to_uniq
                                if next_pos < 0  or next_pos >= len(useful_read):
                                    stopped += 1
                                else:
                                    if not (useful_read[next_pos] in extensions.keys()):
                                        extensions[useful_read[next_pos]] = 0
                                    extensions[useful_read[next_pos]] += 1
                                break
        
                    print ("making desicion...")
                    sum_ext = 0
                    max_ext = ""
                    for extension in extensions:
                        sum_ext += extensions[extension]
                        if max_ext == "" or extensions[extension] > extensions[max_ext]:
                            max_ext = extension
                    if max_ext != "":
                        print (f"extension search stopped {stopped} best {extensions[max_ext]} total {sum_ext}")
                    else:
                        print ("no extension found")
                    if max_ext != "" and extensions[max_ext] * 1.6 >= sum_ext and extensions[max_ext] > 3:
                        print ("going further")
                        progress = True
                        
                        if direction == -1:
                            to_search = to_search + max_ext
                            active_node_list = active_node_list + [max_ext]
                        else:
                            to_search = max_ext + to_search
                            active_node_list = [max_ext] + active_node_list
                        '''if max_ext[1:] in manual_unique:
                            shift_to_uniq = len (active_node_list) -1
                            print (f"Shifting to {shift_to_uniq}")
                            to_search = max_ext'''
                    else:
                        print ("stopped")
                        if max_ext != "" and extensions[max_ext] * 1.6 < sum_ext and extensions[max_ext] > 5:
                            print ("WARN, non-unique continuation")
                            print(extensions)
                        print ("".join(to_search))
                        break
                    if direction == -1:
                        addon[direction] = active_node_list[reuse:]
                    else:
                        addon[direction] = active_node_list[:len(active_node_list) - reuse]
            res_nodes = []
            for j in range (0, len(nodes)):
                if j == i:
                    adleft = 0
                    if (-1 in addon.keys()):
                        for k in range (0, len (addon[-1])):
                            res_nodes.append(addon[-1][k])
                            adleft += nodelens[addon[-1][k]]
                            adleft -= link_len[res_nodes[len(res_nodes) - 2]][res_nodes[len(res_nodes) - 1]]                    
                res_nodes.append(nodes[j])
                if j == i:
                    adright = 0
                    if (1 in addon.keys()):
                        for k in range (0, len (addon[1])):
                            res_nodes.append(addon[1][k])
                            adright += nodelens[addon[1][k]]
                            if k != 0:
                                adright -= link_len[res_nodes[len(res_nodes) - 2]][res_nodes[len(res_nodes) - 1]]                    
                if j == i + 1:
                    if (1 in addon.keys()) and len(addon[1]) > 0:
                        adright -= link_len[res_nodes[len(res_nodes) - 2]][res_nodes[len(res_nodes) - 1]]                                            
            new_seq = "".join(res_nodes)
            print("\t".join([arr[0], new_seq, arr[2]]))
            print (f"Added compressed length {adleft} {adright}")
            intersect = ""
            if (1 in addon.keys()) and (-1 in addon.keys()):
                for chl in addon[-1]:
                    for chr in addon[1]:
                        if chl == chr:
                            intersect += chl      
                        break            
            if intersect == "":
                to_print = "empty intersection " 
            else:
                to_print = "WARN NONEMPTY intersection " + intersect
            print (to_print)
            sys.stdout.flush()
                            
                    
