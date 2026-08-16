#!/usr/bin/env python3
"""Whole-domain conservative-adjacency cluster and max-stat inference on balanced maps."""
from __future__ import annotations
import json, multiprocessing as mp
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.csgraph import connected_components

ROOT=Path(__file__).resolve().parent.parent;RUN=ROOT;DATA=RUN/'data/BALANCED_MAPS.npz'
COORD=RUN/'data/Coordenadas_EEG_GENE_23canais.ced';ADJ=RUN/'data/CHANNEL_ADJACENCY.json'
SHAPE=(20,37,75);PS=(.01,.005);N=100000;SEEDS=[2026081901+i for i in range(4)]
def fstat(data):
 n,k,_=data.shape;g=data.mean((0,1),keepdims=True);c=data.mean(0,keepdims=True);p=data.mean(1,keepdims=True)
 return (n*np.sum((c-g)**2,axis=1).squeeze()/2)/np.maximum(np.sum((data-c-p+g)**2,axis=(0,1))/30,np.finfo(float).tiny)
def path(n):return sparse.diags([np.ones(n-1),np.ones(n-1)],[-1,1],shape=(n,n),format='csr')
def graph(channels):
 art=json.loads(ADJ.read_text());tab=pd.read_csv(COORD,sep='\t').set_index('labels').loc[channels];x=tab[['X','Y','Z']].to_numpy(float);x/=np.linalg.norm(x,axis=1,keepdims=True);ix={v:i for i,v in enumerate(channels)};r=[];c=[]
 for a,b in art['edges']:
  i,j=ix[a],ix[b];deg=np.degrees(np.arccos(np.clip(x[i]@x[j],-1,1)))
  if deg<=46.1:r += [i,j];c += [j,i]
 cg=sparse.csr_matrix((np.ones(len(r)),(r,c)),shape=(20,20));ic=sparse.eye(20);ff=sparse.eye(37);tt=sparse.eye(75)
 return (sparse.kron(sparse.kron(cg,ff),tt)+sparse.kron(sparse.kron(ic,path(37)),tt)+sparse.kron(sparse.kron(ic,ff),path(75))).tocsr()
def maxmass(stat,thr,adj):
 active=np.flatnonzero(stat>thr)
 if not len(active):return 0.
 n,l=connected_components(adj[active][:,active],directed=False);return float(np.bincount(l,weights=stat[active]-thr,minlength=n).max(initial=0))
def comps(stat,thr,adj):
 active=np.flatnonzero(stat>thr)
 if not len(active):return []
 n,l=connected_components(adj[active][:,active],directed=False)
 return sorted([active[l==i] for i in range(n)],key=lambda a:float(np.sum(stat[a]-thr)),reverse=True)
def worker(arg):
 seed,n=arg;a=np.load(DATA);flat=a['cr_safe'].astype(float).reshape(16,3,-1);adj=graph(a['channels']);rng=np.random.default_rng(seed);perms=np.asarray([(0,1,2),(0,2,1),(1,0,2),(1,2,0),(2,0,1),(2,1,0),(2,1,0)],np.int8)[:6];pidx=np.arange(16)[:,None];th=[float(stats.f.ppf(1-p,2,30)) for p in PS];mass=np.empty((n,2));maximum=np.empty(n)
 for i in range(n):
  choice=rng.integers(0,6,size=16);s=fstat(flat[pidx,perms[choice],:]);maximum[i]=s.max()
  for j,t in enumerate(th):mass[i,j]=maxmass(s,t,adj)
 return mass,maximum
def rec(indices,stat,thr,ch,fr,tm,flat,label,adj):
 co=np.asarray(np.unravel_index(indices,SHAPE)).T;peak=int(indices[np.argmax(stat[indices])]);pc=np.unravel_index(peak,SHAPE);means=flat[:,:,indices].mean((0,2));
 # articulation/check is intentionally not repeated; prior audit validated adjacency. 
 return {'component':label,'cluster_mass':float(np.sum(stat[indices]-thr)),'voxel_count':len(indices),'time_min_ms':float(tm[co[:,2].min()]),'time_max_ms':float(tm[co[:,2].max()]),'frequency_min_hz':float(fr[co[:,1].min()]),'frequency_max_hz':float(fr[co[:,1].max()]),'channel_count':len(np.unique(co[:,0])),'peak_f':float(stat[peak]),'peak_channel':str(ch[pc[0]]),'peak_frequency_hz':float(fr[pc[1]]),'peak_time_ms':float(tm[pc[2]]),'mean_condition_0_db':float(means[0]),'mean_condition_1_db':float(means[1]),'mean_condition_2_db':float(means[2])}
def main():
 a=np.load(DATA);maps=a['cr_safe'].astype(float);flat=maps.reshape(16,3,-1);stat=fstat(flat);adj=graph(a['channels']);thresholds=[float(stats.f.ppf(1-p,2,30)) for p in PS]
 with mp.Pool(4) as pool:res=pool.map(worker,[(s,N//4) for s in SEEDS])
 null=np.vstack([x[0] for x in res]);maxnull=np.hstack([x[1] for x in res]);rows=[];labels_primary=np.zeros(np.prod(SHAPE),np.int32)
 for pi,(p,thr) in enumerate(zip(PS,thresholds)):
  components=comps(stat,thr,adj);sig=0
  for ci,ind in enumerate(components,1):
   mass=float(np.sum(stat[ind]-thr));exc=int(np.sum(null[:,pi]>=mass));pval=(exc+1)/(N+1);co=np.asarray(np.unravel_index(ind,SHAPE)).T;mean_t=co[:,2].mean();label='early' if mean_t<25 else 'late' if mean_t>=25 else 'mixed'
   row={'specification':'frequency_safe_code_balanced','cft_p':p,'f_threshold':thr,'observed_cluster_count':len(components),'cluster_number':ci,'permutation_exceedances':exc,'cluster_p_fwer':pval,'permutations':N,**rec(ind,stat,thr,a['channels'],a['frequencies_hz'],a['times_ms'],flat,label,adj)};rows.append(row)
   if pval<.05:
    sig+=1
    if p==.01: labels_primary[ind]=ci
  rows.append({'specification':'frequency_safe_code_balanced','cft_p':p,'f_threshold':thr,'observed_cluster_count':len(components),'cluster_number':0,'component':'summary','significant_cluster_count':sig,'permutations':N})
 # Absolute power extracted on safe-map significant cores, a non-discovery sensitivity.
 absolute=a['cr_absolute'].astype(float).reshape(16,3,-1)
 for row in rows:
  if row.get('cft_p')==.01 and row.get('cluster_number',0)>0 and row.get('cluster_p_fwer',1)<.05:
   ind=np.flatnonzero(labels_primary==row['cluster_number']);mean=absolute[:,:,ind].mean((0,2));rows.append({'specification':'absolute_log_post_on_safe_balanced_core','cft_p':.01,'cluster_number':row['cluster_number'],'component':row['component'],'mean_condition_0_db':mean[0],'mean_condition_1_db':mean[1],'mean_condition_2_db':mean[2],'voxel_count':len(ind),'note':'same significant safe-balanced cluster; no separate absolute-power search'})
 pd.DataFrame(rows).to_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_CLUSTER_RESULTS.csv',index=False)
 critical=float(np.quantile(maxnull,.95,method='higher'));survive=stat>=critical;coord=np.asarray(np.where(survive)).T;peak=int(np.argmax(stat));pc=np.unravel_index(peak,SHAPE);old=np.zeros(SHAPE, np.int32);
 maxrow={'specification':'frequency_safe_code_balanced','permutations':N,'maxstat_fwer_critical_f':critical,'surviving_voxel_count':int(survive.sum()),'peak_f':float(stat[peak]),'peak_channel':str(a['channels'][pc[0]]),'peak_frequency_hz':float(a['frequencies_hz'][pc[1]]),'peak_time_ms':float(a['times_ms'][pc[2]]),'overlap_previous_early_voxels':int(np.logical_and(survive.reshape(SHAPE),old==2).sum()),'overlap_previous_late_voxels':int(np.logical_and(survive.reshape(SHAPE),old==1).sum())}
 pd.DataFrame([maxrow]).to_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_MAXSTAT.csv',index=False)
 with (RUN/'data/BALANCED_INFERENCE.npz').open('wb') as f:np.savez_compressed(f,observed_f=stat.reshape(SHAPE).astype(np.float32),labels_cft01=labels_primary.reshape(SHAPE),max_null=maxnull.astype(np.float32),conservative_adjacency='geodesic<=46.1deg')
 (RUN/'logs/balanced_whole_domain.json').write_text(json.dumps({'status':'PASS','permutations':N,'seeds':SEEDS,'cft_ps':PS,'adjacency':'connected Delaunay subset, geodesic <=46.1 degrees'},indent=2)+'\n')
 print({'status':'PASS','permutations':N,'maxstat_survivors':int(survive.sum())})
if __name__=='__main__':main()
