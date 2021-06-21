from __future__ import print_function,division
import re,codecs,random,csv,json
from alignment import Needleman, Hirschberg
from nltk.tokenize import sent_tokenize
from query_data_from_mongodb import query_abstract_text
import pandas as pd

def extract_function_section_from_uniprot(protein_id_list,uniprot_kb_file):
    #uniprot_kb_file='./Q9NZC2.txt'
    uniprot_function_dic={}
    with codecs.open(uniprot_kb_file,'r','utf8') as kbf:
        count=0
        previous_line_is_ac=False
        for line in kbf:
            line=line.strip()

            if len(line)==0:
                continue
            if line.startswith('AC   '):
                #if there are many AC's for one protein
                #only extract the first one and ignore the others
                if previous_line_is_ac:
                    continue
                previous_line_is_ac=True
                id_list=line[3:].split(';')
                id_list=[i.strip() for i in id_list if i!='']

                for idi in id_list:
                    if idi in protein_id_list:
                        uniprot_function_dic[idi]=''

                #print(id_list)
            elif line.startswith('CC   -!- FUNCTION:'):
                previous_line_is_ac=False
                #sometimes, there will be a section in {}
                for_function_section=True


                for idi in id_list:
                    if idi in protein_id_list:
                        uniprot_function_dic[idi]+=line[len('CC   -!- FUNCTION:'):]

            elif line.startswith('CC       '):
                previous_line_is_ac=False
                if for_function_section:

                    for idi in id_list:
                        if idi in protein_id_list:
                            uniprot_function_dic[idi]+=' '
                            uniprot_function_dic[idi]+=line[len('CC       '):]

            else:
                previous_line_is_ac=False
                for_function_section=False

    json_file='./uniprot_protein_function.json'
    with open(json_file,'w') as outfile:
        json.dump(uniprot_function_dic,outfile)

    #sometimes, there will be a section in {}, delete that
    for fi in uniprot_function_dic.keys():
        function_section_text=uniprot_function_dic[fi]
        if function_section_text.find('}.'):
            left_parathesis_index=function_section_text.rfind('{')
            function_section_text=function_section_text[:left_parathesis_index]
        uniprot_function_dic[fi]=function_section_text


    #divide the function section into sentences by pmid for each protein
    protein_pubmed_dic={}
    for ki in uniprot_function_dic.keys():
        function_section_text=uniprot_function_dic[ki]
        pubmed_dic={}
        while len(function_section_text)>1:
            #print(function_section_text)
            pubmed_index=function_section_text.find('PubMed:')
            by_similarity_index=function_section_text.find('(By similarity')

            #if nothing is found, just break
            if pubmed_index==-1:
                break

            #deal with the cases when there are multiple PubMed for one sentence
            if pubmed_index==0:
                next_pmid=''
                digit_index=pubmed_index+7
                while True:
                    next_letter=function_section_text[digit_index]
                    if next_letter.isdigit():
                        next_pmid+=next_letter
                        digit_index+=1
                    else:
                        break
                if next_pmid not in pubmed_dic:
                    pubmed_dic[next_pmid]=[current_sentence]
                else:
                    pubmed_dic[next_pmid].append(current_sentence)
                function_section_text=function_section_text[digit_index+2:]
                function_section_text=function_section_text.strip()
                continue

            current_sentence=function_section_text[:pubmed_index-1]
            current_sentence=current_sentence.strip()

            if by_similarity_index<pubmed_index and by_similarity_index>0:
                #skip the by similarity cases
                function_section_text=function_section_text[by_similarity_index+len('(By similarity).'):]
                function_section_text=function_section_text.strip()
            else:

                next_pmid=''
                digit_index=pubmed_index+7
                while True:
                    next_letter=function_section_text[digit_index]
                    if next_letter.isdigit():
                        next_pmid+=next_letter
                        digit_index+=1
                    else:
                        break
                if next_pmid not in pubmed_dic:
                    pubmed_dic[next_pmid]=[current_sentence]
                else:
                    pubmed_dic[next_pmid].append(current_sentence)
                function_section_text=function_section_text[digit_index+2:]
                function_section_text=function_section_text.strip()
        protein_pubmed_dic[ki]=pubmed_dic
    #print the unfound protein ids
    for pi in protein_id_list:
        if pi not in protein_pubmed_dic:
            print("Protein not found: ",pi)

    return protein_pubmed_dic
'''
for pi in pubmed_dic.keys():
    sent_list=query_abstract_text(pi)
    print('PMID:',pi)
    print('Abstract:',' '.join(sent_list))
    for ai in pubmed_dic[pi]:
        print('Summary sentence:',ai)
        score_list=[]
        seqa=list(ai)
        for seqb in sent_list:
            seqb=list(seqb)
            n = Needleman()
            a,b = n.align(seqa, seqb)

            score = n.score(a, b)
            #print(score)
            score_list.append(score)
        print('Score for this sentence:',score_list)
'''

def merge_sentence_based_on_pmid(protein_pubmed_dic_json):
    with open(protein_pubmed_dic_json) as jfile:
        protein_pubmed_dic=json.load(jfile)
    pubmed_dic={}
    for proteini in protein_pubmed_dic.keys():
        for pmidi in protein_pubmed_dic[proteini]:
            if pmidi not in pubmed_dic:
                pubmed_dic[pmidi]=[]
            for senti in protein_pubmed_dic[proteini][pmidi]:
                if senti not in pubmed_dic[pmidi] and senti.find('ECO:')<0:
                    pubmed_dic[pmidi].append(senti)
    json_file='./pubmed_dic.json'
    with open(json_file,'w') as outfile:
        json.dump(pubmed_dic,outfile)

if __name__ == "__main__":
    protein_list_file='protein_list.txt'
    proteinList = pd.read_csv(protein_list_file,header=None).iloc[:,0].tolist()

    uniprot_kb_file='uniprot_sprot.dat'
    #protein_pubmed_dic=extract_function_section_from_uniprot(proteinList,uniprot_kb_file)
    #print(protein_pubmed_dic)
    json_file='./uniprot_protein_pubmed_dic.json'
    #with open(json_file,'w') as outfile:
    #    json.dump(protein_pubmed_dic,outfile)
    merge_sentence_based_on_pmid(json_file)


