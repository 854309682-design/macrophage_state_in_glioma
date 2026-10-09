"""
Step 8 — Druggability of TAM-IS candidate genes (DGIdb drugs + Open Targets tractability).

DGIdb is queried for known drug-gene interactions; Open Targets (GraphQL v4) provides
tractability buckets (small-molecule / antibody). Builds the "intervenability" argument
that replaces the (dropped) MR layer.

Outputs: results/tables/step08_druggability.tsv

Run: env/.venv/bin/python scripts/08_causal_druggability/step08_druggability.py
"""
import json
import os
import time
import urllib.request

import pandas as pd

# --- paths: portable — override with GLIOMA_PROJ / GLIOMA_DATA; defaults derive from this file's location ---
_HERE = os.path.dirname(os.path.abspath(__file__))          # scripts/<step>/
ROOT = os.environ.get("GLIOMA_PROJ") or os.path.dirname(os.path.dirname(_HERE))
TBL = os.path.join(ROOT, "results/tables")
OT = "https://api.platform.opentargets.org/api/v4/graphql"
DGIDB = "https://dgidb.org/api/graphql"


def post(url, payload):
    req = urllib.request.Request(url, data=json.dumps(payload).encode(),
                                 headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


TAMIS_CANDIDATES = ["SPP1", "TREM2", "APOE", "APOC1", "C1QA", "C1QB", "LGALS3",
                    "CTSL", "GPNMB", "CCL3", "CD163", "LPL", "TYROBP", "CSF1R", "ITGAM"]

if __name__ == "__main__":
    # ---- DGIdb: one batched query ----
    genes_json = json.dumps(TAMIS_CANDIDATES)
    dg = post(DGIDB, {"query": "{ genes(names:%s){ nodes{ name interactions{ drug{name} } } } }" % genes_json})
    drugs_by_gene = {}
    for node in dg.get("data", {}).get("genes", {}).get("nodes", []) or []:
        names = sorted({i["drug"]["name"] for i in node.get("interactions", []) if i.get("drug")})
        drugs_by_gene[node["name"].upper()] = names

    # ---- Open Targets tractability ----
    def ot_tract(gene):
        try:
            s = post(OT, {"query": 'query($q:String!){search(queryString:$q,entityNames:["target"],page:{index:0,size:1}){hits{id name}}}', "variables": {"q": gene}})
            hits = s["data"]["search"]["hits"]
            if not hits or hits[0]["name"].upper() != gene.upper():
                return ""
            tid = hits[0]["id"]
            t = post(OT, {"query": 'query($id:String!){target(ensemblId:$id){tractability{label modality value}}}', "variables": {"id": tid}})
            tr = t["data"]["target"]["tractability"] or []
            on = [f"{x['modality']}:{x['label']}" for x in tr if x.get("value")]
            return "; ".join(on)
        except Exception:
            return ""

    rows = []
    for g in TAMIS_CANDIDATES:
        d = drugs_by_gene.get(g, [])
        tract = ot_tract(g)
        rows.append({"gene": g, "n_drugs_dgidb": len(d),
                     "top_drugs": ", ".join(d[:6]), "tractability_OT": tract})
        print(f"{g}: DGIdb drugs={len(d)} | OT tractability='{tract[:70]}'", flush=True)
        time.sleep(0.2)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TBL, "step08_druggability.tsv"), sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 70)
    print("\n=== druggability summary ===", flush=True)
    print(df.to_string(index=False), flush=True)
    print("DONE", flush=True)
