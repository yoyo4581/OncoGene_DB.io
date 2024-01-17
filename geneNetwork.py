from subprocess import Popen, PIPE
import time
import streamlit as st
import networkx as nx 
import matplotlib.pyplot as plt

mycursor = st.connection('mysql', type='sql')

def GeneQuery(Gene_sub, Gene_pool, limit, count, iter):
    time.sleep(0.1)
    mycursor.execute("select gene1, gene2 FROM ggLink WHERE (gene1=(%s) OR gene2=(%s)) AND NOT linkTypes='text' AND minResCount<50 ORDER BY pairCount DESC LIMIT %s",(Gene_sub, Gene_sub, int(limit)))
    myresult = mycursor.fetchall()
    print(myresult)
    if count==1:
        return myresult
    for tup in myresult: # go through the first result
        if tup[0] != Gene_sub: #if the target gene is not the same gene as the subject
            retest = GeneQuery(tup[0], Gene_pool,limit, 1, iter) #call the same function again, this time with that target gene
            for tup2 in retest: #iterate through result
                if tup2[0] and tup2[1] in Gene_pool and tup2 not in myresult: #look if any target or source is a gene from our gene pool, if it is add the combination. prevent non unique additions
                    myresult.append(tup2)
        if tup[1] != Gene_sub:
            retest = GeneQuery(tup[1], Gene_pool,limit, 1, iter) #call the same function again, this time with that target gene
            for tup2 in retest: #iterate through result
                if tup2[0] and tup2[1] in Gene_pool and tup2 not in myresult: #look if any target or source is a gene from our gene pool, if it is add the combination. and statement to prevent recall of same subject
                    myresult.append(tup2)
    print(myresult)
    return myresult

def Graphit(node_list):
    G = nx.MultiDiGraph()
    for nodes in node_list:
        G.add_edge(nodes[0], nodes[1])
    return G.reverse()
