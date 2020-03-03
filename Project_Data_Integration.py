#!/usr/bin/env python
# coding: utf-8

# <font size="12">Information Integration Bioinformatics Project</font>
# 

# ## Settings

# In[1]:


UNIPROT_URL = 'https://www.uniprot.org/uniprot/?query=cancer%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22&format=txt&force=true&sort=score'
NCBI_URL = 'ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia'
SAVE_XLSX = False
DOWNLOAD_FILE = True

PSQL_USER = 'adrie'
PSQL_PASSWORD = '123'
PSQL_HOST = '127.0.0.1'
PSQL_PORT = '5432'
PSQL_DATABASE = 'adrie'


# In[2]:


import numpy as np
import pandas as pd
import re
import requests
import urllib.request as request
import gzip
import psycopg2
from psycopg2 import Error


# # Uniprot Part

# In[3]:


if DOWNLOAD_FILE:
    url = UNIPROT_URL
    r = requests.get(url)  
    with open("uniprot.txt",'wb') as f: 
        f.write(r.content)


# In[4]:


file = open("uniprot.txt", 'r') 
lines = file.readlines() 


# In[5]:


uniprot = pd.DataFrame(columns=['ID','GN_name','AC','Lenght','Rec_Name','MIM','KW'])
uni_synonyms = pd.DataFrame(columns=['ID','GN_name','GN_synonym'])
uni_altName = pd.DataFrame(columns=['ID','Alt_name','Short_Alt_name'])


# In[6]:


def make_tables(ID, GN_Name, AC, Lenght, Rec_Name, MIM, KW, GN_Synonyms, alts, alt_short, uniprot, uni_synonyms, uni_altName):
    #uniprot
    uniprot = uniprot.append(pd.DataFrame([[ID, GN_Name, AC, Lenght, Rec_Name, MIM, KW]], 
                    columns=['ID','GN_name','AC','Lenght','Rec_Name','MIM','KW']), ignore_index=True)
    #synonyms
    for syn in GN_Synonyms:
        uni_synonyms = uni_synonyms.append(pd.DataFrame([[ID, GN_Name, syn]], 
                    columns=['ID','GN_name','GN_synonym']), ignore_index=True)
    #altName
    for alt in alts:
        short = ''
        if alt in alt_short:
            short = alt_short[alt]
        uni_altName = uni_altName.append(pd.DataFrame([[ID, alt, short]], 
                    columns=['ID','Alt_name','Short_Alt_name']), ignore_index=True)
        
    return uniprot, uni_synonyms, uni_altName


# In[7]:


def ID_handler(line):
    sp = line.split()
    ID = sp[1]
    Lenght = sp[3]
    return ID, Lenght

def AC_handler(line):
    AC = ''
    sp = line.split()
    for i in range(1,len(sp)):
        AC=AC+sp[i]
    return AC

def DE_handler(line):
    result = re.search('=(.*);', line)
    return result.group(1)

def GN_handler(line, scenario):  # 1.Name=/2.Synonyme=/3.Both
    
    GN_Name = ''
    GN_Synonyms = []
    
    if scenario == 1:
        line = line.split(';')[0]
        sp = line[line.find('Name=')+5:].split()
        for s in sp:
            if ('{' not in s):
                if('ECO' not in s):
                    GN_Name = s
        return GN_Name
        
        
    if scenario == 2:
        sys = re.search('Synonyms=(.*)', line)
        sp = sys.group(1).split()
        
        for s in sp:
            if ('{' not in s) or ('ECO' not in s):
                s = s.replace(',' , '')
                s = s.replace(';', '')
                GN_Synonyms.append(s)
        return GN_Synonyms
        
        
    if scenario == 3:
        GN_Name = line[line.find('Name=')+5:line.find(';')]
        GN_Name = GN_Name.split()[0]
        
        sys_line = line[line.find('Synonyms=')+9:]
        sys_line = sys_line[:sys_line.find(';')]
        sp = sys_line.split()
        for s in sp:
            s = s.replace(',' , '')
            GN_Synonyms.append(s)

    return GN_Name, GN_Synonyms

def MIM_handler(line): # if gene if MIM
    sp = line.split()[2]
    MIM = sp.split(';')[0]
    return MIM

def KW_handler(line): #results should be concated
    KW = line.replace('KW   ' , '')
    return KW


# In[8]:


for line in lines: 
    
    if line.startswith('ID'):
        KW = ''
        ID, Lenght = ID_handler(line)
    
    if line.startswith('AC'):
        AC = AC_handler(line)
        
    if line.startswith('DE'):
        if 'RecName' in line:
            alts = []
            alt_short = dict()
            Rec_Name = DE_handler(line)
        if 'AltName' in line:
            Alt_name = DE_handler(line)
            alts.append(Alt_name)
        if 'Short' in line:
            alt_short[Alt_name] = DE_handler(line)
            
    if line.startswith('GN'):
        scenario = 0 
        GN_Synonyms = []
        if (' Name' in line) and ('Synonyms' not in line):
            scenario = 1
            GN_Name = GN_handler(line, scenario)
        if ('Synonyms' in line) and ('Name' not in line):
            scenario = 2
            GN_Synonyms = GN_handler(line, scenario)
        if (' Name' in line) and ('Synonyms' in line):
            scenario = 3
            GN_Name, GN_Synonyms = GN_handler(line, scenario)

    if line.startswith('DR'):
        if 'MIM' in line:
            if 'gene' in line:
                MIM = MIM_handler(line)

    if line.startswith('KW'):
        KW = KW + KW_handler(line)
            
    if line.startswith('SQ'): #object finished ==> make tables
        uniprot, uni_synonyms, uni_altName = make_tables(ID, GN_Name, AC, Lenght, Rec_Name, MIM, KW, GN_Synonyms, 
                                                 alts, alt_short, uniprot, uni_synonyms, uni_altName)
         


# In[9]:


uniprot.head()


# In[10]:


uni_synonyms.head()


# In[11]:


uni_altName.head()


# In[12]:


if SAVE_XLSX:
    uniprot.to_excel('uniprot_main.xlsx')
    uni_synonyms.to_excel('uniprot_synonyms.xlsx')
    uni_altName.to_excel('uniprot_altName.xlsx')


# # NCBI Part

# In[13]:


human_genome_file = './Homo_sapiens.gene_info'

if DOWNLOAD_FILE:
    
    baseURL = NCBI_URL
    filename = "Homo_sapiens.gene_info.gz"

    response = request.urlopen("/".join([baseURL, filename, ]))
    with open(human_genome_file, "wb") as f:
        f.write(gzip.decompress(response.read()))
        
gene_info_file = human_genome_file
human_genome_info = pd.read_csv(gene_info_file, delimiter='\t')


# Not cleaned table from dataset

# In[14]:


human_genome_info.head()


# Needed columns for the project

# In[15]:


columns = ['GeneID', 'Symbol', 'chromosome', 'map_location', 'description', 'type_of_gene', 'Symbol_from_nomenclature_authority', 'Full_name_from_nomenclature_authority']
select_df = human_genome_info[columns]
select_df.head()


# As we need only MIM value, we extract it from dbXrefs column

# In[16]:


xref_rows = list()

for i, series in human_genome_info.iterrows():
    gene_id = series.GeneID
    symbol = series.Symbol
    dbXrefs = series.dbXrefs
    if pd.notnull(dbXrefs):
        a = dbXrefs.split('|')
        if a[0].startswith('MIM:'):
            b,c=a[0].split(':', 1)
            xref_rows.append((gene_id, symbol, c))
        else:
            xref_rows.append((gene_id, symbol, ' '))
            
xref_df = pd.DataFrame(xref_rows, columns=['GeneID','Symbol', 'MIM'])
xref_df.head()


# Combine two tables to get needed information from table

# In[17]:


nbi_main=pd.merge(select_df, xref_df, on='Symbol')
nbi_main = nbi_main.rename(columns={"Symbol": "GN_Name"})
nbi_main.head()


# In[18]:


synonym_rows = list()
for i, series in human_genome_info.iterrows():
    gene_id = series.GeneID
    symbol = series.Symbol
    synonyms = series.Synonyms
    if pd.notnull(synonyms):
        for synonym in synonyms.split('|'):
            synonym_rows.append((gene_id, symbol,synonym))
nbi_synonym = pd.DataFrame(synonym_rows, columns=['GeneID', 'GN_Name', 'GN_Synonym'])
nbi_synonym.head()


# In[19]:


if SAVE_XLSX:
    nbi_main.to_excel('NBI_main.xlsx')
    xref_df.to_excel('NBI_mim.xlsx')
    nbi_synonym.to_excel('NBI_synonym.xlsx')


# # SQL Part

# In[36]:


def create_tables(user,password,host,port,database):
    try:
        connection = psycopg2.connect(user = user,
                                      password = password,
                                      host = host,
                                      port = port,
                                      database = database)
        cursor = connection.cursor()
       
        query = '''DROP TABLE IF EXISTS Uniprot_Name;
                   DROP TABLE IF EXISTS Uniprot_Synonym;
                   DROP TABLE IF EXISTS Uniprot;
                   DROP TABLE IF EXISTS NBI_Synonym;
                   DROP TABLE IF EXISTS NBI;'''
        cursor.execute(query)
        connection.commit()
        
            
        query = '''
        CREATE TABLE Uniprot(
            ID TEXT,
            GN_name TEXT NOT NULL,
            AC TEXT NOT NULL,
            Length INT NOT NULL,
            Rec_Name TEXT NOT NULL,
            MIM TEXT,
            KW  TEXT NOT NULL,
            PRIMARY KEY (ID)
        );
        CREATE TABLE Uniprot_Name(
            ID TEXT,
            Alt_name TEXT  NOT NULL,
            Short_Alt_name  TEXT,
            FOREIGN KEY (ID) REFERENCES Uniprot(ID)
        );
        CREATE TABLE Uniprot_Synonym(
            ID TEXT,
            GN_name TEXT NOT NULL,
            GN_synonym TEXT NOT NULL,
            FOREIGN KEY (ID) REFERENCES Uniprot(ID)
        );
        CREATE TABLE NBI( 
            ID TEXT,
            GN_name TEXT NOT NULL,
            Chromosone TEXT NOT NULL,
            Map_Location TEXT NOT NULL,
            Description TEXT NOT NULL,
            GN_Type TEXT NOT NULL,   
            NA_Symbol TEXT NOT NULL,  
            NA_FullName TEXT NOT NULL,
            GN_ID INT NOT NULL,
            MIM TEXT,
            PRIMARY KEY (ID)            
        );
        CREATE TABLE NBI_Synonym(
            ID TEXT,
            GN_Name TEXT NOT NULL,
            GN_Synonym TEXT NOT NULL,
            FOREIGN KEY (ID) REFERENCES NBI(ID)
        );'''

        cursor.execute(query)
        connection.commit()
        print("Table created successfully in PostgreSQL ")

    except (Exception, psycopg2.DatabaseError) as error :
        print ("Error while creating PostgreSQL table", error)
    finally:
        #closing database connection.
            if(connection):
                cursor.close()
                connection.close()


# In[29]:


def column_table(dataframe):
    liste = []
    for c in dataframe:
        liste.append(c)
    return liste


def update(user,password,host,port,database,table,dataframe):
    try:
        connection = psycopg2.connect(user = user,
                                      password = password,
                                      host = host,
                                      port = port,
                                      database = database)
        cursor = connection.cursor()
        nb_column = len(dataframe.columns) - 1
        column = column_table(dataframe)
        for index, row in dataframe.iterrows():
            query = "INSERT INTO "+table+" VALUES ("+"%s,"*nb_column+"%s) ON CONFLICT DO NOTHING"
            cursor = connection.cursor()
            data = ()
            for c in column:
                data +=(row[c],)
            cursor.execute(query, data)
            connection.commit()
        print ("Records inserted successfully into "+table+" table")
    except (Exception, psycopg2.DatabaseError) as error :
        print ("Error while creating PostgreSQL table", error)
    finally:
         if(connection):
                cursor.close()
                connection.close()


# In[30]:


def select(user,password,host,port,database,query):
    try:
        connection = psycopg2.connect(user = user,
                                      password = password,
                                      host = host,
                                      port = port,
                                      database = database)
        cursor = connection.cursor()
        cursor.execute(query)
        connection.commit()
        records = cursor.fetchall()
        print("Total rows: ", len(records))
        print("Printing each row")
        for row in records:
            result = ""
            for x in range(len(row)):
                result += str(row[x])+"\t "
            print(result)
        cursor.close()
    except (Exception, psycopg2.DatabaseError) as error :
        print ("Error while creating PostgreSQL table", error)
    finally:
         if(connection):
                cursor.close()
                connection.close()
                print("PostgreSQL connection is closed")


# In[37]:


create_tables(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE) 
update(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, "uniprot", uniprot)
update(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, "uniprot_synonym", uni_synonyms)
update(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, "uniprot_name", uni_altName)
update(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, "NBI", nbi_main)
update(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, "NBI_Synonym", nbi_synonym)


# ## SQL queries

# ### Genes with one common name in NCBI and UniProtKB but with different other names in NCBI and UniProtKB

# In[24]:


query = "SELECT DISTINCT U.GN_name FROM Uniprot_Synonym U, NBI_Synonym N WHERE U.GN_NAME=N.GN_NAME AND U.GN_Synonym!=N.GN_Synonym"
select(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, query)


# Total rows:  4717
# Printing each row
# SMAD4	 
# WWOX	 
# ADGRF4	 
# FAM83D	 
# MLLT1	 
# PARD6A	 
# CPSF4	 
# CLIP1	 
# MSL1	 
# ZNF142	 
# AASDH 	 
# .	 
# .	 
# .

# ### Genes with one common name in NCBI and UniProtKB but where the MIM terms differs between NCBI and UniProtKB

# In[25]:


query = "SELECT DISTINCT U.GN_name FROM Uniprot U, NBI N WHERE U.GN_NAME=N.GN_NAME AND U.MIM!=N.MIM"
select(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, query)


# Total rows:  908
# Printing each row
# CPXCR1	 
# WIZ	 
# MAPK1IP1L	 
# SH3BP5L	 
# PRPF38B	 
# INO80D	 
# RAB3IL1	 
# TTC29	 
# RPS6KL1	 
# LENG1	
# .	 
# .	 
# .

# ### All the alternative names for the gene 'UBE2S_HUMAN'

# In[26]:


query = "SELECT DISTINCT N.GN_synonym FROM Uniprot U,NBI_Synonym N WHERE U.ID='UBE2S_HUMAN' AND N.GN_Name=U.GN_Name UNION SELECT DISTINCT US.GN_synonym FROM Uniprot U, Uniprot_Synonym US WHERE U.ID='UBE2S_HUMAN' AND US.GN_Name=U.GN_Name"
select(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, query)


# ### Genes with common MIM name in NCBI and UniProtKB 

# In[27]:


query = "SELECT DISTINCT U.GN_name, U.MIM FROM Uniprot U, NBI N WHERE U.MIM=N.MIM UNION SELECT DISTINCT N.GN_name, N.MIM FROM Uniprot U, NBI N WHERE U.MIM=N.MIM"
select(PSQL_USER, PSQL_PASSWORD, PSQL_HOST, PSQL_PORT, PSQL_DATABASE, query)


# Total rows:  7738
# Printing each row
# AADAC	 600338	 
# AAK1	 616405	 
# AAMP	 603488	 
# AARS	 601065	 
# AARS1	 601065	 
# AARSD1	 613212	 
# AASDH	 614365	 
# AASDHPPT	 607756	 
# AATF	 608463	 
# AATK	 605276	
# .	 
# .	 
# .
