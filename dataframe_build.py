import streamlit as st

# Python Program for demonstrating the  
# PyMongo Cursor to Pandas DataFrame 
  
   
# Importing required modules 
 
from pandas import DataFrame 
import pymongo
import pandas as pd
  
   
# Connecting to MongoDB server 
# client = MongoClient('host_name', 
# 'port_number') 
def init_connection():
    return pymongo.MongoClient(st.secrets["mongo"].uri)

client = init_connection()
  
# Connecting to the database named 
# GFG 
cancer_db = client.Cancer_db 
   
# Accessing the collection named 
# gfg_collection 
labels = cancer_db.Labels 
query = cancer_db.Query_data
genes = cancer_db.Gene_data
# Now creating a Cursor instance 
# using find() function 
cursor1 = labels.find()
cursor2 = query.find()
cursor3 = genes.find()
print('Type of cursor:',type(cursor1)) 
  
# Converting cursor to the list of  
# dictionaries 
list_cur1 = list(cursor1)
list_cur2 = list(cursor2)
list_cur3 = list(cursor3)

for item in list_cur2:
    item['GEO']=item['value']['GEO']
    if item['value']['PubMedIds']:
        item['PubMedIds'] = item['value']['PubMedIds']
    elif not item['value']['PubMedIds']:
        item['PubMedIds'] = None
    if item['value']['title']:
        item['title'] = item['value']['title']
    elif not item['value']['title']:
        item['title'] = None

  
# Converting to the DataFrame 
df1 = DataFrame(list_cur1)
df2 = DataFrame(list_cur2) #multiple dictionaries in this list -> df
df3 = DataFrame(list_cur3)


df2 = df2.rename(columns={"_id":"_id2"})
df2 = df2.drop('value', axis=1)

df12 = pd.merge(df1, df2, left_on='_id', right_on='link', how='left').drop('_id', axis=1)

df12 = df12.drop('link', axis=1)

df123 = pd.merge(df12, df3, left_on='_id2', right_on='link', how='left').drop('_id2', axis=1)

df123 = df123.drop('link', axis=1)

print('Type of df:',type(df1), type(df2), type(df123)) 
  
# Printing the df to console 
print() 
print(df1.head())

print() 
print(df2.head())

print()
print(df12.head())

print()
print(df123.head())
