from __future__ import print_function,division
import re,codecs,random,csv,json
from alignment import Needleman, Hirschberg
from nltk.tokenize import sent_tokenize
from query_data_from_mongodb import query_abstract_text

#uniprot_kb_file='./Glygen/OGER/uniprot_sprot_human.dat'
uniprot_kb_file='./Q9NZC2.txt'
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
            #id_list=id_list[:1]
            if 'Q9NZC2' in id_list or 'O43914' in id_list:
                print(id_list)
                for idi in id_list:
                    uniprot_function_dic[idi]=''

            #print(id_list)
        elif line.startswith('CC   -!- FUNCTION:'):
            #sometimes, there will be a section in {}
            for_function_section=True
            
            if 'Q9NZC2' in id_list or 'O43914' in id_list:
                for idi in id_list:
                    uniprot_function_dic[idi]+=line[len('CC   -!- FUNCTION:'):]

        elif line.startswith('CC       '):
            if for_function_section:
                if 'Q9NZC2' in id_list or 'O43914' in id_list:
                    for idi in id_list:
                        uniprot_function_dic[idi]+=' '
                        uniprot_function_dic[idi]+=line[len('CC       '):]

        else:
            for_function_section=False

#print(uniprot_function_dic['Q9NZC2'])
#sometimes, there will be a section in {}, delete that
for fi in uniprot_function_dic.keys():
    function_section_text=uniprot_function_dic[fi]
    if function_section_text.endswith('}.'):
        left_parathesis_index=function_section_text.rfind('{')
        function_section_text=function_section_text[:left_parathesis_index]
    uniprot_function_dic[fi]=function_section_text

#print(uniprot_function_dic['Q9NZC2'])

#divide the function section into sentences by pmid
function_section_text=uniprot_function_dic['Q9NZC2']
pubmed_dic={}
while len(function_section_text)>1:
    #print(function_section_text)
    pubmed_index=function_section_text.find('PubMed:')
    by_similarity_index=function_section_text.find('(By similarity')
    
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

print(pubmed_dic)
#json_file='./Glygen/OGER/uniprot_human_sprot.json'
#with open(json_file,'w') as outfile:
#    json.dump(uniprot_kb_dic,outfile)

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
