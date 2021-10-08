import pandas as pd
import re

pairs_file = '/home/mvazque2/.rbbt/var/jobs/ExTRI/pairs/Default.tsv'
output_file = "evidence_counts_ExTRI_pairs.tsv"

tsv = pd.read_csv(pairs_file, sep='\t', header=1, dtype='string')
columns = tsv.columns

db_pmid_c = {}
db_sign_c = {}

for c in columns:
    m = re.search(r'\[(.*)\] PMID', c)
    if m and m.group(1):
        db = m.group(1)
        db_pmid_c[db] = c

    m = re.search(r'\[(.*)\] Effect', c)
    if m and m.group(1):
        db = m.group(1)
        if not db in db_sign_c:
            db_sign_c[db] = c

    m = re.search(r'\[(.*)\] Regulation', c)
    if m and m.group(1):
        db = m.group(1)
        db_sign_c[db] = c

    m = re.search(r'\[(.*)\] Activation/Repression', c)
    if m and m.group(1):
        db = m.group(1)
        db_sign_c[db] = c

    m = re.search(r'\[(.*)\] Sign', c)
    if m and m.group(1):
        db = m.group(1)
        db_sign_c[db] = c

dbs = db_pmid_c.keys()

pair_evidence = {}
for i, row in tsv.iterrows():
    pair = row["#TF:TG"]
    evidence = [[], [], []]
    for db in dbs:

        pmid = row[db_pmid_c[db]]
        if row.isna()[db_pmid_c[db]]:
            pmid = ""
        pmid = pmid.split("|")

        if db in db_sign_c:
            sign = row[db_sign_c[db]]
            if row.isna()[db_sign_c[db]]:
                sign = ""
            sign = sign.split("|")
        else:
            sign = ["Unknown"] * len(pmid)

        for i in zip(pmid, sign):
            pl = i[0]

            ps = pl.split(";")

            s = i[1]

            if s == "":
                s = "Unknown"
            elif s == "Activation":
                s = "UP"
            elif s == "positive":
                s = "UP"
            elif s == "Repression":
                s = "DOWN"
            elif s == "negative":
                s = "UP"
            elif s == "unknown":
                s = "Unknown"
            elif s == "Unknown":
                s = "Unknown"
            elif s == "UP":
                s = "UP"
            elif s == "DOWN":
                s = "DOWN"

            for p in ps:
                if p != "":
                    if s == "UP":
                        evidence[0].append(p)
                    elif s == "DOWN":
                        evidence[1].append(p)
                    elif s == "Unknown":
                        evidence[2].append(p)
                    else:
                        print("ERROR: ", s)

            pair_evidence[pair] = list(map(lambda x: len(list(set(x))), evidence))

output = pd.DataFrame.from_dict(pair_evidence, orient='index', columns=["UP", "DOWN", "Unknown"])
output = output.rename_axis("#Pair")

output.to_csv(output_file, sep="\t")
