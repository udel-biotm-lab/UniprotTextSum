from __future__ import division,print_function
from collections import OrderedDict
from pymongo import MongoClient
from bson.son import SON
from bson.codec_options import CodecOptions
#import pandas as pd
import re,csv,os,json

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
    sentence_list=[]
    if raw_doc is not None:
        #print(raw_doc)
        sentence = raw_doc["sentence"]

        for senInfo in sentence:
            #senText_original = title_abstract[senInfo["charStart"]:senInfo["charEnd"]]
            sentence_list.append(raw_doc['text'][int(senInfo["charStart"]):int(senInfo["charEnd"])+1])

    return sentence_list, raw_doc

def generate_csv_labeling(pubmed_dic,uniprot_id_to_gene_name,output_file):

    csv_title=['score','gene','pmid','summary_sent_index','abstract_sent_index','comments','summary_sent','abstract_sent',]
    #csv_row=['']*len(csv_title)
    pmid_list=[pi for pi in pubmed_dic.keys()]

    row_count=0
    with open(output_file, 'w') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerow(csv_title)
        for pmid in pmid_list:
            sentence_list, doc_info=query_abstract_text(pmid)
            if doc_info is not None:
                sentence = doc_info["sentence"]
                for si in range(len(pubmed_dic[str(pmid)])):

                    for eachSen in sentence:
                        if "index" not in eachSen:
                            index = 0
                        else:
                            index = eachSen["index"]

                        # SenToDisplay[index] = abstractText[eachSen["charStart"]:eachSen["charEnd"]]
                        sentenceText = doc_info['text'][int(eachSen["charStart"]):int(eachSen["charEnd"])+1]

                        csv_row=['']*len(csv_title)
                        csv_row[csv_title.index('gene')]='=HYPERLINK(\"https://www.uniprot.org/uniprot/'+pubmed_dic[str(pmid)][si][1]+'\"'+',\"'+str(uniprot_id_to_gene_name[pubmed_dic[str(pmid)][si][1]])+'\")'  \
                            if pubmed_dic[str(pmid)][si][1] in uniprot_id_to_gene_name else ''
                        csv_row[csv_title.index('pmid')]='=HYPERLINK(\"http://biotm2.cis.udel.edu:8009/pmid/'+str(pmid)+'\"'+',\"'+str(pmid)+'\")'
                        csv_row[csv_title.index('summary_sent_index')]=str(si+1)
                        csv_row[csv_title.index('abstract_sent_index')]=str(index)
                        csv_row[csv_title.index('summary_sent')]=pubmed_dic[str(pmid)][si][0]
                        csv_row[csv_title.index('abstract_sent')]=sentenceText
                        csv_row=[ri.encode('utf-8') for ri in csv_row ]
                        spamwriter.writerow(csv_row)
                        row_count+=1
            #if row_count>10000:
            #    break


if __name__ == "__main__":
    pubmed_dic_json='pubmed_dic.json'
    output_file='for_labeling.csv'
    with open(pubmed_dic_json) as jfile:
        pubmed_dic=json.load(jfile)
    uniprot_id_to_gene_name_json_file='./uniprot_id_to_gene_name.json'
    with open(uniprot_id_to_gene_name_json_file) as jfile:
        uniprot_id_to_gene_name=json.load(jfile)
    generate_csv_labeling(pubmed_dic,uniprot_id_to_gene_name,output_file)

