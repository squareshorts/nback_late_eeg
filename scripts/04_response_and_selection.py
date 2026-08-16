#!/usr/bin/env python3
"""Correctness-selection, false-alarm timing and response-status analyses."""
from __future__ import annotations
import itertools,json
from pathlib import Path
import numpy as np,pandas as pd
from scipy import signal,special,stats
import statsmodels.formula.api as smf
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parent.parent;RUN=ROOT;P7=ROOT/'submission_pipeline/de_novo_analysis/phase7_foundation_lock';P75=ROOT/'submission_pipeline/de_novo_analysis/phase7_5_contract_redteam'
MAN='';TECH='';COORD=RUN/'data/Coordenadas_EEG_GENE_23canais.ced'
FS=600.;WORK=np.arange(-600,1200);POST=600+np.arange(0,900,12);FREQ=np.arange(4.,41.);CHANNELS=np.array(['F7','T3','T5','Fp1','F3','C3','P3','O1','F8','T4','T6','Fp2','F4','C4','P4','O2','Fz','Cz','Pz','Oz']);CODES=np.array([f'C{i}' for i in range(1,11)]);PIDS=np.array([10,16,18,19,22,24,27,28,30,31,33,34,36,38,39,40])
def first(path):
 with path.open(encoding='ascii') as f:
  for x in f:
   if x.strip():return x.strip()
def load(path):
 h=list(CHANNELS[:16])+['Fz','Cz','Pz','Oz','A1','A2','FOTO'];return np.loadtxt(path,delimiter='\t',skiprows=int(first(path).split('\t')==h))
def g(c,m=4):
 out=np.zeros_like(c,float)
 for n in range(1,51):out+=((2*n+1)/(n**m*(n+1)**m))*special.eval_legendre(n,c)
 return out/(4*np.pi)
def interp_weights(x,target):
 target=np.asarray(target,int);good=np.asarray([i for i in range(20) if i not in set(target)],int);gg=g(np.clip(x[good]@x[good].T,-1,1));a=np.block([[gg+1e-5*np.eye(len(good)),np.ones((len(good),1))],[np.ones((1,len(good))),np.zeros((1,1))]]);gt=g(np.clip(x[target]@x[good].T,-1,1));return good,(np.c_[gt,np.ones(len(target))]@np.linalg.inv(a))[:,:-1]
def wave(f):
 cyc=3+5*(f-4)/36;s=cyc/(2*np.pi*f);half=int(np.ceil(4*s*FS));t=np.arange(-half,half+1)/FS;k=np.exp(2j*np.pi*f*t)*np.exp(-t*t/(2*s*s));return k/np.sqrt(np.sum(np.abs(k)**2)),half
def exact(v):
 signs=np.asarray(list(itertools.product((-1.,1.),repeat=16)),np.float32);return float(np.mean(np.abs(signs@v/16)>=abs(v.mean())-1e-12))
def main():
 if (RUN/'data/BALANCED_LOSO_ALL_NT_TRIAL_FEATURES.csv').exists():
  print('Bypassing 04 raw generation')
  features=pd.read_csv(RUN/'data/BALANCED_LOSO_ALL_NT_TRIAL_FEATURES.csv')
 else:
  (RUN/'data').mkdir(exist_ok=True);m=pd.read_csv(MAN);q=pd.read_csv(TECH);d=m.merge(q[['participant_id','trial_index_within_session','flagged_channels','alternative_technical_pass']],on=['participant_id','trial_index_within_session']);nt=d[(~d.target_status_ALVO)&d.alternative_technical_pass].copy();masks=np.load(RUN/'data/BALANCED_LOSO_MASKS.npz')['masks'];tab=pd.read_csv(COORD,sep='\t').set_index('labels').loc[CHANNELS];xyz=tab[['X','Y','Z']].to_numpy(float);xyz/=np.linalg.norm(xyz,axis=1,keepdims=True);sos=signal.butter(4,[.5,45],btype='bandpass',fs=FS,output='sos');rows=[];overlap=[]
  for pi,pid in enumerate(PIDS):
   tr=nt[nt.participant_id.eq(pid)].sort_values('trial_index_within_session').reset_index(drop=True);raw=load(ROOT/tr.raw_file_path.iloc[0]);sample=tr.foto_sample.to_numpy(int);ep=np.empty((len(tr),20,len(WORK)),np.float32)
   for ch in range(20):ep[:,ch]=signal.sosfiltfilt(sos,raw[:,ch],padtype='odd',padlen=27)[sample[:,None]+WORK[None,:]]
   cache={}
   for i,names in enumerate(tr.flagged_channels.fillna('')):
    if names:
     target=tuple(sorted(np.flatnonzero(np.isin(CHANNELS,names.split('|'))).tolist()));
     if target not in cache:cache[target]=interp_weights(xyz,target)
     good,w=cache[target];ep[i,np.asarray(target)]=w@ep[i,good]
   ep-=ep.mean(axis=1,keepdims=True);feat=np.zeros((len(tr),2));voxel_times=[]
   for ci in range(2):voxel_times.append(np.where(masks[pi,ci])[2])
   for fi,f in enumerate(FREQ):
    k,half=wave(float(f));safeidx=np.arange(300,600)[WORK[np.arange(300,600)]<=-half];coef=signal.fftconvolve(ep,k[None,None,:],mode='same',axes=-1);power=np.abs(coef)**2;safe=power[...,safeidx].mean(axis=-1)
    for ci in range(2):
     co=np.asarray(np.where(masks[pi,ci])).T;pick=np.flatnonzero(co[:,1]==fi)
     if len(pick):
      ch=co[pick,0];tm=co[pick,2];post=np.take(power,POST[tm],axis=2)[:,ch,np.arange(len(ch))] if False else None
      # Avoid coupled advanced indexing: select voxel columns one at a time, then accumulate log ratio.
      for j,(cc,tt) in enumerate(zip(ch,tm)):
       feat[:,ci]+=10*np.log10(power[:,cc,POST[tt]]/safe[:,cc])
   feat/=np.array([masks[pi,0].sum(),masks[pi,1].sum()])
   for i,r in tr.iterrows():
    for ci,name in enumerate(('early','late')):rows.append({'participant_id':pid,'load':r.load,'stimulus_code':str(r.stimulus_code).upper(),'behavior_class':r.behavior_class,'response_status':int(r.response_present),'response_latency_ms':r.response_latency_ms,'cluster':name,'trial_ersp_db':feat[i,ci]})
    if r.behavior_class=='NONTARGET_FALSE_ALARM' and np.isfinite(r.response_latency_ms):
     for ci,name in enumerate(('early','late')):
      t=voxel_times[ci]*20.;rt=float(r.response_latency_ms);overlap.append({'participant_id':pid,'load':r.load,'cluster':name,'rt_ms':rt,'mask_voxel_time_n':len(t),'before_response_fraction':float(np.mean(t<rt-200)),'within_plusminus_200_fraction':float(np.mean(np.abs(t-rt)<=200)),'after_response_fraction':float(np.mean(t>rt+200))})
   print(f'response features participant {pid:03d}',flush=True)
  features=pd.DataFrame(rows);features.to_csv(RUN/'data/BALANCED_LOSO_ALL_NT_TRIAL_FEATURES.csv',index=False);ov=pd.DataFrame(overlap);summary=ov.groupby('cluster').agg(false_alarm_trial_n=('rt_ms','size'),rt_median_ms=('rt_ms','median'),rt_q25_ms=('rt_ms',lambda x:np.quantile(x,.25)),rt_q75_ms=('rt_ms',lambda x:np.quantile(x,.75)),before_response_fraction=('before_response_fraction','mean'),within_plusminus_200_fraction=('within_plusminus_200_fraction','mean'),after_response_fraction=('after_response_fraction','mean')).reset_index();summary.to_csv(RUN/'FALSE_ALARM_RT_OVERLAP.csv',index=False);ov.to_csv(RUN/'data/FALSE_ALARM_RT_OVERLAP_TRIALS.csv',index=False)
 modelrows=[]
 for cluster in ('early','late'):
  w=features[features.cluster.eq(cluster)].copy();w['response_status']=w.response_status.astype(int)
  try:
   fit=smf.gee('trial_ersp_db ~ C(load) * response_status',groups='participant_id',data=w,cov_struct=sm.cov_struct.Exchangeable(),family=sm.families.Gaussian()).fit()
   for term in fit.params.index:
    ci=fit.conf_int().loc[term];modelrows.append({'cluster':cluster,'term':term,'estimate_db':fit.params[term],'ci95_low':ci.iloc[0],'ci95_high':ci.iloc[1],'p_value':fit.pvalues[term],'model':'GEE Gaussian exchangeable; participant clusters','trial_n':len(w),'false_alarm_n':int(w.response_status.sum()),'stability':'sparse false alarms; interaction exploratory'})
  except Exception as e:modelrows.append({'cluster':cluster,'term':'MODEL_FAILED','model':str(e),'trial_n':len(w)})
 pd.DataFrame(modelrows).to_csv(RUN/'RESPONSE_STATUS_MODELS.csv',index=False)
 # 500 exact code-balanced, equal-count ALL-NT sensitivity on independent trial features.
 rng=np.random.default_rng(20260820);sens=[]
 for cluster in ('early','late'):
  values=np.empty((500,16,3))
  for pi,pid in enumerate(PIDS):
   w=features[(features.participant_id==pid)&(features.cluster==cluster)&features.stimulus_code.isin(CODES)]
   cells={(lo,c):w[(w.load==lo)&(w.stimulus_code==c)].trial_ersp_db.to_numpy() for lo in range(3) for c in CODES};mins={c:min(len(cells[lo,c]) for lo in range(3)) for c in CODES}
   for rep in range(500):
    for lo in range(3):values[rep,pi,lo]=np.mean(np.concatenate([rng.choice(cells[lo,c],mins[c],replace=False) for c in CODES]))
  step=((values[:,:,1]+values[:,:,2])/2-values[:,:,0]).mean(1);diff=(values[:,:,2]-values[:,:,1]).mean(1)
  for kind,v in [('step_later_minus_0',step),('2_minus_1',diff)]:
   direction = 1 if cluster=='early' and kind.startswith('step') else -1 if cluster=='late' and kind.startswith('step') else np.sign(np.median(v))
   sens.append({'cluster':cluster,'contrast':kind,'repetitions':500,'median_db':float(np.median(v)),'interval_2_5':float(np.quantile(v,.025)),'interval_97_5':float(np.quantile(v,.975)),'percent_expected_direction':float(100*np.mean(np.sign(v)==direction)),'estimator':'trial-level safe ERSP within held-out balanced LOSO masks; equal code and count sampling'})
 pd.DataFrame(sens).to_csv(RUN/'ALL_NT_BALANCED_TRIAL_SENSITIVITY.csv',index=False)
if __name__=='__main__':main()