name_search, name_search_rev, ids = ncbi_gene_search_by_name(cgenes, 'Homo Sapiens')
name_gene_annot = search_gene_annot(ids)

for i in name_gene_annot:
    name_search_rev[i].extend(name_gene_annot[i])

name_search_rev_rev = dict(zip([x[0] for x in name_search_rev.values() if x !=['-']], name_search_rev.values()))
search_names = list(name_search_rev_rev.keys())