#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import division, print_function
import time, re, sys, os, ast
import pandas as pd
import pymongo
import urllib
import json
from pymongo import MongoClient
import random, operator
from collections import OrderedDict
from bson.son import SON
from bson.codec_options import CodecOptions
import json
from xml.etree import ElementTree as ET
from flask import Flask, jsonify, request, render_template, send_from_directory, url_for, Response
from flask import session
from flask import send_file
from nltk.tokenize import word_tokenize

input_uniprot_pubmed_json = './pubmed_dic.json'
with open(input_uniprot_pubmed_json) as jfile:
    pubmed_dic = json.load(jfile)

uniprot_id_to_gene_name_json_file='./uniprot_id_to_gene_name.json'
with open(uniprot_id_to_gene_name_json_file) as jfile:
    uniprot_id_to_gene_name_dic = json.load(jfile)
# --- create database instances---
# Environment variables
mongodb_host = os.environ.get("MONGODB_HOST", "localhost")  # change to biotm2.cis.udel.edu before dockerizing
mongodb_port = os.environ.get("MONGODB_PORT", "27017")
db_name = os.environ.get("DBNAME", "medline_current")  # change database name for your own dbName
textCollectionName = os.environ.get("COLLECTION_TEXT", "text")

# Database URI
MONGODB_URI = 'mongodb://' + mongodb_host + ':' + mongodb_port + '/'

# Database object
client = MongoClient(MONGODB_URI)
opts = CodecOptions(document_class=SON)

# Database
dbName = client[db_name]  # glygen

# Collection
textCollection = dbName[textCollectionName].with_options(codec_options=opts)

# ------------ Static Files------------#
# fixedFilePath = "/data/Applications/debarati/Glygen"
# sys.path.append(fixedFilePath)


# --- create Flask instance---
app = Flask(__name__)
app.config["APPLICATION_ROOT"] = os.environ.get("APPLICATION_ROOT", "labeling/")


# --Route functions
@app.route('/')
def index():

    pmid_list = [ki for ki in pubmed_dic.keys()]

    return render_template('index.html', listOfPmids=pmid_list)


@app.route('/pmid/<pmid>', methods=['POST', 'GET'])
def uniprot(pmid):
    SenToDisplay = OrderedDict()
    SumSenToDisplay = OrderedDict()
    abstractDoc = textCollection.find_one({"docId": str(pmid)})
    if abstractDoc:

        if "sentence" in abstractDoc and "text" in abstractDoc:
            abstractText = abstractDoc["text"]
            senInfo = abstractDoc["sentence"]
            SenToDisplay = OrderedDict()

            for eachSen in senInfo:
                if "index" not in eachSen:
                    index = 0
                else:
                    index = eachSen["index"]

                # SenToDisplay[index] = abstractText[eachSen["charStart"]:eachSen["charEnd"]]
                sentenceText = abstractText[eachSen["charStart"]:eachSen["charEnd"]+1]

                SenToDisplay[index] = sentenceText

        for si in range(len(pubmed_dic[str(pmid)])):
            SumSenToDisplay[si+1]= '<a href=\"https://www.uniprot.org/uniprot/'+pubmed_dic[str(pmid)][si][1]\
            +'#names_and_taxonomy\">'+uniprot_id_to_gene_name_dic[pubmed_dic[str(pmid)][si][1]]+'</a>'+pubmed_dic[str(pmid)][si][0]
        return render_template('showAbstract.html', text=SenToDisplay, pmid=pmid, summary_text=SumSenToDisplay)
    else:
        return render_template('exitPage.html')

if __name__ == "__main__":
    host = sys.argv[1]
    # host = '0.0.0.0'
    port = sys.argv[2]
    is_debug = os.environ.get("IS_DEBUG", "false")
    debugMode = False
    if is_debug.lower() == "true":
        debugMode = True
    if not port:
        port = '4500'
    try:
        app.secret_key = ".."
        # app.secret_key = "1341dea76d6f"
        SESSION_PERMANENT = False
        '''
            from werkzeug.serving import run_simple
            from werkzeug.wsgi import DispatcherMiddleware
            app.config['DEBUG'] = debugMode
            # # Load a dummy app at the root URL to give 404 errors.
            # # Serve app at APPLICATION_ROOT for localhost development.
            application = DispatcherMiddleware(Flask('dummy_app'), {
                app.config['APPLICATION_ROOT']: app,
            })
            run_simple(host, int(port), application, use_reloader=True)
        '''
        app.run(host=host, port=port, debug=debugMode)
    except KeyboardInterrupt:
        print('Interrupted')
        try:
            sys.exit(0)
        except SystemExit:
            os._exit(0)
