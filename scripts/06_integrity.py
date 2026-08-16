#!/usr/bin/env python3
"""Validate required stopping-analysis deliverables and create a hash manifest."""
from pathlib import Path
import csv,hashlib,json
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parent.parent;RUN=ROOT
required='''FINAL_ROBUSTNESS_EXECUTIVE_SUMMARY.md
COMMON_STIMULUS_AVAILABILITY.csv
DETERMINISTIC_STIMULUS_BALANCED_CLUSTER_RESULTS.csv
DETERMINISTIC_STIMULUS_BALANCED_MAXSTAT.csv
STIMULUS_BALANCED_LOSO_RESULTS.csv
ORIGINAL_VS_BALANCED_EFFECTS.csv
NONTARGET_SELECTION_COUNTS.csv
ALL_NONTARGET_EARLY_RESULTS.csv
ALL_NONTARGET_LATE_RESULTS.csv
FALSE_ALARM_RT_OVERLAP.csv
RESPONSE_STATUS_MODELS.csv
HARD_CASE_ANALYSIS.csv'''.splitlines()
missing=[x for x in required if not (RUN/x).is_file()]
if missing:raise RuntimeError(missing)
cluster=pd.read_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_CLUSTER_RESULTS.csv')
if not np.isfinite(cluster.f_threshold.dropna()).all() or not set(cluster.cft_p.dropna())=={.01,.005}:raise RuntimeError('Invalid CFT threshold table')
if int(cluster.permutations.dropna().max())!=100000:raise RuntimeError('Cluster permutations not 100,000')
maxstat=pd.read_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_MAXSTAT.csv')
if int(maxstat.permutations.iloc[0])!=100000 or int(maxstat.surviving_voxel_count.iloc[0])!=12:raise RuntimeError('Maxstat gate')
figures=list((RUN/'figures').glob('*.pdf'))
if len(figures)!=7:raise RuntimeError(f'Expected 7 figures, got {len(figures)}')
manifest=[]
for p in sorted(x for x in RUN.rglob('*') if x.is_file() and x.name!='HASH_MANIFEST.csv' and '__pycache__' not in x.parts):manifest.append({'relative_path':str(p.relative_to(RUN)),'size_bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
with (RUN/'HASH_MANIFEST.csv').open('w',newline='',encoding='utf-8') as f:
 w=csv.DictWriter(f,fieldnames=['relative_path','size_bytes','sha256']);w.writeheader();w.writerows(manifest)
out={'status':'PASS','required_deliverables':len(required),'figures':len(figures),'hashed_files':len(manifest),'permutations':100000,'maxstat_surviving_voxels':12}
(RUN/'logs/integrity.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
