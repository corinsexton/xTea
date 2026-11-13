import sys,glob,gzip,statistics as st
for arm in ("base","v2"):
  for d in sorted(glob.glob(arm+"/*/")):
    fs=[f for f in glob.glob(d+"*/*.vcf")]
    if not fs: continue
    op=gzip.open if fs[0].endswith(".gz") else open
    v={"0/1":[], "1/1":[]}; w={"0/1":[], "1/1":[]}
    for l in op(fs[0],"rt"):
      if l[0]=="#": continue
      t=l.rstrip("\n").split("\t"); k=t[8].split(":"); s=dict(zip(k,t[9].split(":")))
      g=s.get("GT")
      if g in v and s.get("VAF",".") not in (".",""):
        v[g].append(float(s["VAF"]))
        ci=s.get("VAF_CI",".").split(",")
        if len(ci)==2 and "." not in ci: w[g].append(float(ci[1])-float(ci[0]))
    print(arm,d.split("/")[1],*["%s n=%d med=%.3f ciw=%.3f"%(g,len(v[g]),st.median(v[g]),st.median(w[g]) if w[g] else -1) for g in v])
