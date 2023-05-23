
from Bio import Entrez
import numpy as np
import re
import pandas as pd
import os
import urllib.request

Entrez.email = "suhininaev@mail.ru"

df_ids = pd.read_csv('../only_ids.txt', header=None)

for id in df_ids[0]:
    #get summary
    temphandle = Entrez.read(Entrez.esummary(db="assembly", id=id, retmode="xml"))
    url = temphandle['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_GenBank']
    label = os.path.basename(url)
    #get the fasta link - change this to get other formats
    link = os.path.join(url,label+'_genomic.gbff.gz')
    urllib.request.urlretrieve(link, f'{label}.gbff.gz') 
