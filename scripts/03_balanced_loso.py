#!/usr/bin/env python3
"""LOSO masks reconstructed only from deterministic common-code balanced maps."""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import sparse,stats
from scipy.sparse.csgraph import connected_components
import importlib.util

ROOT=Path(__file__).resolve().parent.parent;RUN=ROOT;DATA=RUN/'data/BALANCED_MAPS.npz';SHAPE=(20,37,75)
spec=importlib.util.spec_from_file_location('balanced',RUN/'scripts/02_balanced_whole_domain.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
THRESH=float(stats.f.ppf(.99,2,28));PIDS=np.array([10,16,18,19,22,24,27,28,30,31,33,34,36,38,39,40]);CL=('early','late')
def exact(values):
 s=np.asarray(list(itertools.product((-1.,1.),repeat=16)),np.float32);return float(np.mean(np.abs(s@values/16)>=abs(values.mean())-1e-12))
def select(stat,adj,domain):
 active=np.flatnonzero((stat>THRESH)&domain.ravel());n,l=connected_components(adj[active][:,active],directed=False)
 parts=[active[l==i] for i in range(n)];return max(parts,key=lambda x:float(np.sum(stat[x]-THRESH)))
def fstat(data):return mod.fstat(data)
def one_result(frame,population,cluster):
 d=frame[(frame.population==population)&(frame.cluster==cluster)].pivot(index='participant_id',columns='load',values='value').loc[PIDS]
 for contrast,vals in {'step_later_minus_0':(d[1]+d[2])/2-d[0],'2_minus_1':d[2]-d[1]}.items():
  mean=float(vals.mean());half=float(stats.t.ppf(.975,15)*stats.sem(vals));yield {'population':population,'cluster':cluster,'contrast':contrast,'estimate_db':mean,'ci95_low':mean-half,'ci95_high':mean+half,'paired_dz':mean/float(vals.std(ddof=1)),'exact_signflip_p':exact(vals.to_numpy()),'participant_n':16}
def main():
 a=np.load(DATA);maps={'CR':a['cr_safe'].astype(float),'ALL_NT':a['all_safe'].astype(float),'ALL_NT_RESPONSE_CENSORED':a['all_response_censored_safe'].astype(float)};adj=mod.graph(a['channels']);times=a['times_ms'];domains={'early':np.broadcast_to((times<500)[None,None,:],SHAPE),'late':np.broadcast_to((times>=500)[None,None,:],SHAPE)}
 full_masks={}
 full_stat=fstat(maps['CR'].reshape(16,3,-1))
 for name,domain in domains.items():
  selected=select(full_stat,adj,domain);mask=np.zeros(np.prod(SHAPE),bool);mask[selected]=True;full_masks[name]=mask.reshape(SHAPE)
 folds=np.zeros((16,2,*SHAPE),bool);features=[];stab=[]
 for hold,pid in enumerate(PIDS):
  train=np.delete(maps['CR'],hold,axis=0).reshape(15,3,-1);stat=fstat(train)
  for ci,name in enumerate(CL):
   ind=select(stat,adj,domains[name]);mask=np.zeros(np.prod(SHAPE),bool);mask[ind]=True;mask=mask.reshape(SHAPE);folds[hold,ci]=mask
   dice=2*np.logical_and(mask,full_masks[name]).sum()/(mask.sum()+full_masks[name].sum())
   stab.append({'heldout_participant_id':pid,'cluster':name,'voxel_count':int(mask.sum()),'dice_with_full_balanced_mask':float(dice),'threshold_f':THRESH,'selection_rule':'largest CFT=.01 component in fixed nonoverlapping temporal domain from 15 balanced participants'})
   for population,array in maps.items():
    for lo in range(3):features.append({'participant_id':pid,'load':lo,'population':population,'cluster':name,'value':float(array[hold,lo][mask].mean()),'heldout_mask_independent':True})
 rows=[]
 for pop in maps:
  for name in CL:rows.extend(one_result(pd.DataFrame(features),pop,name))
 result=pd.DataFrame(rows)
 result['whole_domain_primary_cft01_support']=result.cluster.eq('late')
 result['interpretation_note']=np.where(result.cluster.eq('late'),
     'balanced LOSO feature result aligned with the FWER-significant balanced primary whole-domain late cluster',
     'conditional LOSO feature result: balanced broad early component did not pass primary whole-domain CFT=.01 FWER inference')
 result.to_csv(RUN/'STIMULUS_BALANCED_LOSO_RESULTS.csv',index=False);pd.DataFrame(stab).to_csv(RUN/'STIMULUS_BALANCED_LOSO_STABILITY.csv',index=False);pd.DataFrame(features).to_csv(RUN/'data/BALANCED_LOSO_FEATURES.csv',index=False)
 with (RUN/'data/BALANCED_LOSO_MASKS.npz').open('wb') as f:np.savez_compressed(f,masks=folds,participant_ids=PIDS,clusters=np.array(CL),times_ms=times,frequencies_hz=a['frequencies_hz'],channels=a['channels'])
 print(result.to_string(index=False))
if __name__=='__main__':main()
