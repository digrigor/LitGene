import requests
import re
import pandas as pd
import unicodedata
import time
import concurrent
import concurrent.futures

pmc_basename_name = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pmc&term="
pubmed_basename_name = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term="
pubmed_basename_uid = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&linkname=gene_pubmed&id="
gene_basename_uid = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id="
gene_basename_name = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term="

# Retrieve a single page and report the url and contents
def load_url(url, timeout):
    response=requests.get(url)
    return response.content

def read_from_csv(path):
    try:
        gs_pd = pd.read_csv(path, header=None)
    except UnicodeDecodeError:
        gs_pd = pd.read_csv(path, encoding="latin1", header=None)

    gs = list(gs_pd[0])

    try: gs_dec = [unicodedata.normalize("NFKD", x).rstrip() for x in gs]
    except TypeError: gs_dec = gs

    return gs_dec


def ncbi_gene_search_by_name(genes, organism, uid=False, annot=True):

    organism = "%20".join(organism.split(" "))

    urlz= []
    url_dic={}

    for g in genes:
        url = gene_basename_name+'('+g+'%5BGene%20Name%5D)%20AND%20'+organism+"%5BOrganism%5D"
        urlz.append(url)
        url_dic[url] = [g]


    with concurrent.futures.ThreadPoolExecutor(max_workers=60) as executor:
        # Start the load operations and mark each future with its URL
        future_to_url = {executor.submit(load_url, url, 60): url for url in urlz}
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            data = future.result()
            try: cor_nam = re.findall("<Term>\S+Gene Name]</Term>", str(data))[0][6:-18]
            except IndexError: cor_nam = "-"
            try:
                fre = re.findall('<Id>\d+</Id>', str(data))[0]
                idz = re.findall('\d+', fre)[0]
            except IndexError: idz = "-"
            url_dic[url].extend([cor_nam, idz])

    url_dic_rev = dict(zip([x[2] for x in url_dic.values() if x !=['-']], url_dic.values()))
    id_list = set(list(url_dic_rev.keys()))
    id_list = [x for x in id_list if x != '-']
    return(url_dic, url_dic_rev, id_list)



def search_gene_annot(uids):
    urlz = []
    uids2 = []
    ext_dic={}

    for i in uids:
        url=gene_basename_uid+str(i)
        urlz.append(gene_basename_uid+str(i))

    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        # Start the load operations and mark each future with its URL
        future_to_url = {executor.submit(load_url, url, 60): url for url in urlz}
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            data = future.result()
            try: gene_symbols = re.findall('<Name>\S+</Name>', str(data))[0][6:-7]
            except IndexError: gene_symbols = "-"
            try: descriptions = re.findall('<Description>.+</Description>', str(data))[0][13:-14]
            except IndexError: descriptions = "-"
            try: summaries = re.findall('<Summary>.+</Summary>', str(data))[0][9:-10]
            except IndexError: summaries = "-"
            uids2 = re.findall("\d+$", url)[0]
            try: ext_dic[uids2].extend([gene_symbols, descriptions, summaries])
            except KeyError: ext_dic[uids2] = [gene_symbols, descriptions, summaries]

    return(ext_dic)

def ncbi_pubmed_name_comention(genes, coterm):
    coterm = "+".join(coterm.split(" "))
    urlz=[]
    pub_dict={}
    for g in genes:
        url = pubmed_basename_name+'('+g+')+AND+('+coterm+')&retmax=10000000'
        urlz.append(url)
        pub_dict[url] = [g]
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        # Start the load operations and mark each future with its URL
        future_to_url = {executor.submit(load_url, url, 60): url for url in urlz}
        for future in concurrent.futures.as_completed(future_to_url):
            url = future_to_url[future]
            data = future.result()
            fre=re.search('<Count>\d+</Count>', str(data)).group()
            counts = re.findall("\d+", fre)
            fre=re.findall('<Id>\d+</Id>', str(data))
            pub_ids=[', '.join(re.findall("\d+",x)) for x in fre]
            result = counts + pub_ids
            pub_dict[url].extend(result)
    return(pub_dict)


dio = ncbi_pubmed_name_comention(search_names, 'Paracetamol')
urlz = []
for g in genes:
    urlz.append(gene_basename_name+'('+g+'%5BGene%20Name%5D)%20AND%20'+organism+"%5BOrganism%5D")


def ncbi_pubmed_uid_comention(genes, coterm):
    coterm = "%20".join(coterm.split(" "))
    counts = []
    pub_ids= []
    for g in genes:
        if g=="-":
            counts.append("-")
            pub_ids.append("-")

        else:
            url=pubmed_basename_uid+g+'&term=('+coterm+')'
            response=requests.get(url)
            a=response.content
            fre=re.findall('<Link>\S{16}\d+', str(a))
            counts.append(len(fre))
            pub_ids.append([','.join(re.findall("\d+",x)) for x in fre])

    return(counts, pub_ids)


a,b,c,d = ncbi_gene_search(["TP53", "TP63", "TP73", "MPOURDELO", "ERBB2", "DEFA1B", "MALAKAS", "TNFAIP6", "CAMP"], "Homo Sapiens", uid=False)


gs = read_from_csv("./germline_genelist_v2.csv")
start = time.time()
a,b,c,d = ncbi_gene_search(gs, "Homo Sapiens", uid=False)
end = time.time()
print(end - start)
start = time.time()
e,f = ncbi_pubmed_uid_comention(b, "Cervical Cancer")
start = time.time()
e,f = ncbi_pubmed_name_comention(a, "Cervical Cancer")
end = time.time()
print(end - start)
Write_them_down("Cristiana5.csv")

def write_them_down(name):
    df = pd.DataFrame()
    df["name"]=a
    df["ncbi_id"]=b
    df["description"]=c
    df["summary"]=d
    df["pubmed_related_articles_counts"]=[int(x[0]) for x in e]
    df["pubmed_related_articles_ids"]=[", ".join(x) for x in f]
    df.to_csv("./"+name)
    return



    df2 = df[df.pubmed_counts != "-"]
    df2["pubmed_counts"]=df2["pubmed_counts"].astype(float)
    df2["pubmed_counts_l"] = df2["pubmed_counts"]+1
    df2["pubmed_counts_l"] = np.log(df2["pubmed_counts_l"])
    ordered_df = df2.sort_values(by='pubmed_counts_l')
    my_range=range(1,len(df2.index)+1)

i=0
for rect in bar:
    height=list(ordered_df["pubmed_counts"])[i]
    plt.text(rect.get_y() + rect.get_width()/2.0, height, '%d' % int(height), ha='center', va='bottom')
    i+=1



# We can use a with statement to ensure threads are cleaned up promptly
nam_id = {}

        nam_id[nam]=id
