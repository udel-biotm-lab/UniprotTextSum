from __future__ import division,print_function
from collections import OrderedDict
from pymongo import MongoClient
from bson.son import SON
from bson.codec_options import CodecOptions
#import pandas as pd
import re,csv,os

def query_abstract_text(pmid):

    mongodb_host = os.environ.get("MONGODB_HOST","0.0.0.0") # change to biotm2.cis.udel.edu before dockerizing
    mongodb_port = os.environ.get("MONGODB_PORT","27017")
    db_name = os.environ.get("DBNAME","medline_current") # change database name for your own dbName
    textCollectionName = os.environ.get("COLLECTION_TEXT","text")

    # Database URI
    MONGODB_URI = 'mongodb://'+mongodb_host+':'+mongodb_port+'/'

    # Database object
    client = MongoClient(MONGODB_URI)
    opts = CodecOptions(document_class=SON)

    # Database
    dbName = client[db_name] # glyco/unicarb

    # Collection
    textCollection = dbName[textCollectionName].with_options(codec_options=opts)

    raw_doc = textCollection.find_one({"docId":str(pmid)})

    sentence = raw_doc["sentence"]
    sentence_list=[]
    for senInfo in sentence:
        #senText_original = title_abstract[senInfo["charStart"]:senInfo["charEnd"]]
        sentence_list.append(raw_doc['text'][int(senInfo["charStart"]):int(senInfo["charEnd"])+1])

    return sentence_list








