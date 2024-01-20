import dataframe_build
import streamlit as st
import pandas as pd
import operator
import geneNetwork
import time
import networkx as nx 
import matplotlib.pyplot as plt
from matplotlib.pyplot import text
import math
from time import sleep
from streamlit_plotly_events import plotly_events
import plotly.graph_objects as go

#page layout

st.set_page_config(layout="wide")
st.title("🐙 OncoGene_db with Graphical Gene Network")
st.subheader("A curated compilation of cancer data from GEO, Literature Web Parsing, and UCSC gene interaction network")


def value_change(category):
    st.session_state.df = st.session_state.df[st.session_state.df[category]==st.session_state[category]]

#builds the dataframe from the database only once
@st.cache_data
def get_data():
    df = dataframe_build.df123
    df = df.loc[:, df.columns!="UID"]
    return df

df = get_data()

#creates a sidebar with dependent multiselect widgets based on available data
with st.sidebar:
    option1 = st.multiselect("Primary Cancer Type", list(set(df.iloc[:, 0])))

    option2 = st.multiselect("Secondary Cancer Type", list(set(df.loc[(df.iloc[:, 0].isin(option1))][df.columns[1]])))

    option3 = st.multiselect("Classification", list(set(df.loc[(df.iloc[:, 0].isin(option1)) & (df.iloc[:, 1].isin(option2))][df.columns[2]])))

    option4 = st.multiselect("Location", list(set(df.loc[(df.iloc[:, 0].isin(option1)) & (df.iloc[:, 1].isin(option2)) & (df.iloc[:, 2].isin(option3))][df.columns[3]])))

    option5 = st.multiselect("Methodology", list(set(df.loc[(df.iloc[:, 0].isin(option1)) & (df.iloc[:, 1].isin(option2)) & (df.iloc[:, 2].isin(option3)) & (df.iloc[:, 3].isin(option4))][df.columns[4]])))

    button_clicked = st.button("Edit dataframe")
    button_reset = st.button('Reset')

#builds the dataframe as a session_state object
if 'df' not in st.session_state or button_reset:
    st.session_state.df = df


#edits data frame
if button_clicked:
    if option1 and option2 and option3 and option4 and option5:
        st.session_state.df = st.session_state.df.loc[
            (st.session_state.df['Primary Cancer Type'].isin(option1)) &
            (st.session_state.df['Secondary Cancer Type'].isin(option2)) &
            (st.session_state.df['Classification'].isin(option3)) &
            (st.session_state.df['Location'].isin(option4)) &
            (st.session_state.df['Method'].isin(option5))]
        button_clicked=False

#converts the dataframe into a downloadable csv file.
def convert_df(df):
    return df.to_csv(index=False).encode('utf-8')

#converts the network into a downloadable csv file.
def convert_Network(list_oftuples):
    net_df = pd.DataFrame(list_oftuples, columns=['target', 'source'])
    return net_df.to_csv(index=False).encode('utf-8')

with st.sidebar:
    csv = convert_df(st.session_state.df)
    button_download = st.download_button('Download Dataframe', csv, "file.csv", "text/csv", key='download-csv')




#user_choice = list(product(*[option1, option2, option3]))

#user_choice = list(*[option1, option2, option3])

#st.write(option1, option2, option3, option4, option5)

#col_part1, col_part2 = st.columns([3,2])

abs_dict = {}
#creates the main data frame and expander widget
#with col_part1:
current_selection = []
with st.container():
    st.dataframe(st.session_state.df.loc[:, ~df.columns.isin(["Primary Cancer Type", "Secondary Cancer Type", "Classification", "Location", "Method"])])

buttons = []

if 'Gene Profile' not in st.session_state:
    st.session_state['Gene Profile']= False

with st.expander('Gene Profile'):

    abstract = st.session_state.df['Abstract'].tolist()
    abstract = [x for x in abstract if x !=[]]

    body = st.session_state.df['Body'].tolist()
    body = [x for x in body if x !=[]]

    keywords = st.session_state.df['Keywords'].tolist()
    keywords = [x for x in keywords if x !=[]]

    abbrev = st.session_state.df['Abbreviations'].tolist()
    abbrev = [x for x in abbrev if x !=[]]

    def plotCall(string):
        st.session_state['current_state'] = string

    generate_genes = st.checkbox("generate Genes", key='generateGenes')
    gene_dict = {}

    if 'generateGenes' not in st.session_state:
        st.session_state.generateGenes = False

    #keeps tract of the frequency at which each gene pops up across contexts
    if generate_genes:
        for abs in abstract:
            if str(abs) !='nan':
                for ab in abs:
                    if str(ab) not in gene_dict:
                        gene_dict[str(ab)] = 1
                    else:
                        gene_dict[str(ab)]+=1
        for bod in body:
            if str(bod) !='nan':
                for item_body in bod:
                    if str(item_body) not in gene_dict:
                        gene_dict[str(item_body)] = 1
                    else:
                        gene_dict[str(item_body)]+=1
        for keyw in keywords:
            if str(keyw) !='nan':
                for key in keyw:
                    if str(key) not in gene_dict:
                        gene_dict[str(key)] = 1
                    else:
                        gene_dict[str(key)]+=1
        for abr in abbrev:
            if str(abr) !='nan':
                for item_abr in abr:
                    if str(item_abr) not in gene_dict:
                        gene_dict[str(item_abr)] = 1
                    else:
                        gene_dict[str(item_abr)]+=1
        abs_dict = dict(sorted(gene_dict.items(), key=operator.itemgetter(1), reverse=True))

        gene_names = [gene for gene in abs_dict.keys()]
        gene_freq = [abs_dict[name] for name in gene_names]


        empty = st.empty()
        with empty.container():
            fig = go.Figure(data=[go.Pie(labels=gene_names, values=gene_freq, hole=0.3)])
            fig.update_traces(hoverinfo='label+value+name', textinfo='none')
            fig.update_layout(uniformtext_minsize=12, margin=dict(t=0, b=0, l=0, r=0),uniformtext_mode='hide')
            selected_points = plotly_events(fig)

        def updatePi(geneList):
            # Reset figure
            FIG = go.Figure()
            indexList = [gene_names.index(i) for i in geneList]
            zero_array = [0]*len(gene_freq)
            zero_array[indexList[-1]] = 0.6  #last one in list will be selected
            for i in indexList:
                if i != indexList[-1]:
                    zero_array[i] = 0.2
    
            Gene_select = geneList[-1]
            st.text(zero_array)
            # Streamlit
            with empty.container():
                st.write("d")
                FIG = go.Figure(data=[go.Pie(labels=gene_names, values=gene_freq, pull=zero_array, hole=0.3)])
                FIG.update_traces(hoverinfo='label+value+name', textinfo='none')
                FIG.update_layout(uniformtext_minsize=12, margin=dict(t=0, b=0, l=0, r=0))
                selected_points = plotly_events(FIG)
                plotCall(Gene_select)
            sleep(1)

        if selected_points:
            point = selected_points[0]['pointNumber']
            if 'selected' not in st.session_state or st.session_state['selected'] != point:
                st.session_state['selected'] = point
                st.write(st.session_state['selected'])
                updatePi([gene_names[point]])

#method to check and display if a gene belongs in the gene pool, or it is new
#also checks if there is gene interaction data for the gene term
def check_result(result):
    covered = []
    in_parse = []
    if result:
        sub_col1, sub_col2 = st.columns([1,1])
        with sub_col1:
            st.subheader("Genes included in literature Parse")
        with sub_col2:
            st.subheader("Genes excluded in literature Parse")
        for tup in result:
            with sub_col1:
                if tup[0] in abs_dict.keys() and tup[0] not in covered:
                    st.markdown(f":red[{str(tup[0])}] \n")
                    covered.append(tup[0])
                    in_parse.append(tup[0])
                if tup[1] in abs_dict.keys() and tup[1] not in covered:
                    st.markdown(f":red[{str(tup[1])}] \n")
                    covered.append(tup[1])
                    in_parse.append(tup[1])
            with sub_col2:
                if tup[0] not in abs_dict.keys() and tup[0] not in covered:
                    st.markdown(str(tup[0]))
                    covered.append(tup[0])
                if tup[1] not in abs_dict.keys() and tup[1] not in covered:
                    st.markdown(str(tup[1]))
                    covered.append(tup[1])
        in_parse.append(st.session_state['current_state'])
        updatePi(in_parse)

    else:
        st.subheader('No Gene data')

#generates graph and information about genes in gene network. Also has a button to download gene network.
#with col_part2:
if 'Network Search' not in st.session_state:
    st.session_state['Network Search']= False
with st.expander('Network Search'):
    if 'current_state' in st.session_state:
        st.header('Subject selected: '+st.session_state['current_state'])
        values = st.slider('Select a number of neighbors to consider in search', 0, 10, 5)
        expansiveText = st.checkbox("Expansive Text", key='expandText')
        if 'expandText' not in st.session_state:
            st.session_state.expandText = False
        generate_graph = st.button('Generate Graph')

        #once the scope has been set and the generate graph button is clicked, then send an sql query
        if generate_graph: 
            st.subheader(f'Top {values} results of genes associated with subject and their associations with genes in gene pool. Pathway database substantiated results only')
            if abs_dict:
                result = geneNetwork.GeneQuery(st.session_state['current_state'], list(abs_dict.keys()), values, 0, 0)
                if result: #if there is a valid query result
                    #generate graph object
                    Graph = geneNetwork.Graphit(result)
                    d = dict(Graph.degree)
                    fig, ax = plt.subplots()
                
                    #try different graph layout if information can be displayed in a simple planar layout 
                    try:
                        pos = nx.planar_layout(Graph)
                    except:
                        st.write("Data not planar, showing spring_layout")
                        pos = nx.spring_layout(Graph)
                        if expansiveText:
                            nx.draw(Graph, pos=pos, with_labels=False, node_size = [(v+math.log(v))*100 for v in d.values()])
                            for node, (x,y) in pos.items():
                                text(x,y, node, fontsize=math.log(d[node]+1)*8,ha='center', va='center')
                        else:
                            nx.draw(Graph, pos=pos, with_labels=False)
                            for node, (x,y) in pos.items():
                                text(x,y, node,ha='center', va='center')
                        st.pyplot(fig)
                        check_result(result)
                    else: #otherwise display information in spring-layout
                        if expansiveText:
                            nx.draw(Graph, pos=pos, with_labels=False, node_size = [(v+math.log(v))*100 for v in d.values()])
                            for node, (x,y) in pos.items():
                                text(x,y, node, fontsize=math.log(d[node]+1)*8, ha='center', va='center')
                        else:
                            nx.draw(Graph, pos=pos, with_labels=False)
                            for node, (x,y) in pos.items():
                                text(x,y, node, ha='center', va='center')
                        st.pyplot(fig)
                        check_result(result)
                    net = convert_Network(result)
                    button_Net = st.download_button('Download Network', net, "Netfile.csv", "text/csv", key='download-csv1')
                else: #if there is no valid gene interaction data display a message for the user
                    st.markdown(f":red[Suspect gene has no interaction data] \n")
