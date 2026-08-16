#!/usr/bin/env python3
"""Construct deterministic common-code balanced maps before feature discovery.

For each participant, condition and recoverable Presentation identity C1--C10,
the code-specific power estimator is computed first.  The condition map is the
unweighted mean of the ten code-specific maps; no random subsampling occurs.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import signal, special

ROOT=Path(__file__).resolve().parent.parent
RUN=ROOT
P7=ROOT/'submission_pipeline/de_novo_analysis/phase7_foundation_lock'
P75=ROOT/'submission_pipeline/de_novo_analysis/phase7_5_contract_redteam'
MANIFEST=''
TECH=''
COORD=RUN/'data/Coordenadas_EEG_GENE_23canais.ced'
FS=600.; WORK=np.arange(-600,1200); POST=600+np.arange(0,900,12); FREQ=np.arange(4.,41.)
CHANNELS=np.array(['F7','T3','T5','Fp1','F3','C3','P3','O1','F8','T4','T6','Fp2','F4','C4','P4','O2','Fz','Cz','Pz','Oz'])
PIDS=np.array([10,16,18,19,22,24,27,28,30,31,33,34,36,38,39,40]); CODES=np.array([f'C{i}' for i in range(1,11)])

def first(path):
    with path.open(encoding='ascii') as f:
        for x in f:
            if x.strip():return x.strip()
    raise RuntimeError(path)
def raw_load(path):
    header=list(CHANNELS[:16])+['Fz','Cz','Pz','Oz','A1','A2','FOTO']; skip=int(first(path).split('\t')==header)
    x=np.loadtxt(path,delimiter='\t',skiprows=skip)
    if x.ndim!=2 or x.shape[1]!=23:raise RuntimeError(f'raw shape {x.shape}')
    return x
def g(c,m=4,terms=50):
    out=np.zeros_like(c,float)
    for n in range(1,terms+1):out+=((2*n+1)/(n**m*(n+1)**m))*special.eval_legendre(n,c)
    return out/(4*np.pi)
def xyz():
    t=pd.read_csv(COORD,sep='\t').set_index('labels').loc[CHANNELS]; x=t[['X','Y','Z']].to_numpy(float)
    return x/np.linalg.norm(x,axis=1,keepdims=True)
def weights(x,target):
    target=np.asarray(target,int);good=np.asarray([i for i in range(20) if i not in set(target)],int)
    gg=g(np.clip(x[good]@x[good].T,-1,1)); a=np.block([[gg+1e-5*np.eye(len(good)),np.ones((len(good),1))],[np.ones((1,len(good))),np.zeros((1,1))]])
    gt=g(np.clip(x[target]@x[good].T,-1,1));return good,(np.c_[gt,np.ones(len(target))]@np.linalg.inv(a))[:,:-1]
def wave(f):
    cycles=3+5*(f-4)/36;sigma=cycles/(2*np.pi*f);half=int(np.ceil(4*sigma*FS));t=np.arange(-half,half+1)/FS
    return np.exp(2j*np.pi*f*t)*np.exp(-t*t/(2*sigma*sigma))/np.sqrt(np.sum(np.exp(-t*t/(sigma*sigma)))),half
def save_npz(path,**data):
    with path.with_suffix('.tmp').open('wb') as f:np.savez_compressed(f,**data)
    path.with_suffix('.tmp').replace(path)

def main():
    (RUN/'data').mkdir(parents=True,exist_ok=True);(RUN/'logs').mkdir(exist_ok=True)
    if (RUN/'data/BALANCED_MAPS.npz').exists():
        print('Bypassing 01 raw generation')
        return
    m=pd.read_csv(MANIFEST); q=pd.read_csv(TECH)
    d=m.merge(q[['participant_id','trial_index_within_session','flagged_channels','alternative_technical_pass']],on=['participant_id','trial_index_within_session'],validate='one_to_one')
    nt=d[(~d.target_status_ALVO)&d.alternative_technical_pass].copy()
    # All maps use shared C identities; P identities are intentionally excluded.
    x=xyz();sos=signal.butter(4,[.5,45],btype='bandpass',fs=FS,output='sos')
    shape=(16,3,20,37,75)
    cr_safe=np.empty(shape,np.float32);cr_abs=np.empty(shape,np.float32);all_safe=np.empty(shape,np.float32);all_abs=np.empty(shape,np.float32);all_censored_safe=np.empty(shape,np.float32)
    avail=[];selection=[]
    for pi,pid in enumerate(PIDS):
        trials=nt[nt.participant_id.eq(pid)].sort_values('trial_index_within_session').reset_index(drop=True)
        raw=raw_load(ROOT/trials.raw_file_path.iloc[0]);samples=trials.foto_sample.to_numpy(int);epochs=np.empty((len(trials),20,len(WORK)),np.float32)
        for ch in range(20):
            filt=signal.sosfiltfilt(sos,raw[:,ch],padtype='odd',padlen=27);epochs[:,ch]=filt[samples[:,None]+WORK[None,:]]
        del raw
        cache={}
        for ti,names in enumerate(trials.flagged_channels.fillna('')):
            if not names:continue
            target=tuple(sorted(np.flatnonzero(np.isin(CHANNELS,names.split('|'))).tolist()))
            if target not in cache:cache[target]=weights(x,target)
            good,w=cache[target];epochs[ti,np.asarray(target)]=w@epochs[ti,good]
        epochs-=epochs.mean(axis=1,keepdims=True)
        code=trials.stimulus_code.str.upper().to_numpy(); load=trials.load.to_numpy(int); iscr=trials.behavior_class.eq('NONTARGET_CORRECT_REJECTION').to_numpy()
        response=trials.response_present.to_numpy(bool);rt=trials.response_latency_ms.to_numpy(float)
        # Sums indexed population(CR/ALL), load, code, channel, frequency, time.
        post=np.zeros((2,3,10,20,37,75),np.float64);base=np.zeros((2,3,10,20,37),np.float64);count=np.zeros((2,3,10),int)
        cens_post=np.zeros((3,10,20,37,75),np.float64);cens_count=np.zeros((3,10,75),int)
        for ci,c in enumerate(CODES):
            for lo in range(3):
                member=(code==c)&(load==lo); selection.append({'participant_id':pid,'load':lo,'stimulus_code':c,'population':'ALL_NT','technical_n':int(member.sum()),'included':bool(member.any())})
                crmember=member&iscr;selection.append({'participant_id':pid,'load':lo,'stimulus_code':c,'population':'CR','technical_n':int(crmember.sum()),'included':bool(crmember.any())})
                count[0,lo,ci]=crmember.sum();count[1,lo,ci]=member.sum()
        for fi,f in enumerate(FREQ):
            kernel,half=wave(float(f));safe_idx=np.arange(300,600)[WORK[np.arange(300,600)]<=-half]
            coef=signal.fftconvolve(epochs,kernel[None,None,:],mode='same',axes=-1);power=np.abs(coef)**2;safe=power[...,safe_idx].mean(axis=-1)
            for ci,c in enumerate(CODES):
                for lo in range(3):
                    member=(code==c)&(load==lo);crmember=member&iscr
                    for pop,idx in enumerate((crmember,member)):
                        if not idx.any():continue
                        post[pop,lo,ci,:,fi]+=np.take(power[idx], POST, axis=2).sum(axis=0)
                        base[pop,lo,ci,:,fi]+=safe[idx].sum(axis=0)
                    # false alarms may contribute only before RT-200; CR trials always contribute.
                    for ti in np.flatnonzero(member):
                        keep=np.ones(75,bool) if not response[ti] else ((POST-600)/FS*1000 < rt[ti]-200)
                        if keep.any():
                            current = cens_post[lo,ci,:,fi,:]
                            current[:,keep] += np.take(power[ti], POST[keep], axis=1)
                            cens_post[lo,ci,:,fi,:] = current
                            cens_count[lo,ci,keep]+=1
        if np.any(count==0):raise RuntimeError(f'common-code cell missing participant {pid}: {np.argwhere(count==0).tolist()}')
        for lo in range(3):
            for pop,(safeout,absout) in enumerate(((cr_safe,cr_abs),(all_safe,all_abs))):
                mp=post[pop,lo]/count[pop,lo,:,None,None,None];mb=base[pop,lo]/count[pop,lo,:,None,None]
                safeout[pi,lo]=np.mean(10*np.log10(mp/mb[...,None]),axis=0).astype(np.float32)
                absout[pi,lo]=np.mean(10*np.log10(mp),axis=0).astype(np.float32)
            mp=cens_post[lo]/np.maximum(cens_count[lo,:,None,None,:],1)
            # Censored estimator retains code-specific safe baselines; no time point is silently imputed.
            mb=base[1,lo]/count[1,lo,:,None,None]
            all_censored_safe[pi,lo]=np.mean(10*np.log10(mp/mb[...,None]),axis=0).astype(np.float32)
        for lo in range(3):avail.append({'participant_id':pid,'load':lo,'common_code_count_K':10,'common_codes':'|'.join(CODES),'missing_common_code_cells':0})
        print(f'balanced maps participant {pid:03d}: {len(trials)} technically valid nontargets',flush=True)
    pd.DataFrame(avail).to_csv(RUN/'COMMON_STIMULUS_AVAILABILITY.csv',index=False)
    pd.DataFrame(selection).to_csv(RUN/'CODE_CELL_COUNTS.csv',index=False)
    # Population counts / correctness selection are reported separately and include P identities as valid nontargets.
    countrows=[]
    for (pid,lo),cell in nt.groupby(['participant_id','load']):
        cr=int(cell.behavior_class.eq('NONTARGET_CORRECT_REJECTION').sum());fa=int(cell.behavior_class.eq('NONTARGET_FALSE_ALARM').sum());alln=len(cell)
        countrows.append({'participant_id':pid,'load':lo,'technical_valid_nontargets':alln,'correct_rejection_nontargets':cr,'false_alarm_nontargets':fa,'other_excluded_nontargets':alln-cr-fa,'percent_excluded_by_correctness':100*(alln-cr)/alln})
    pd.DataFrame(countrows).to_csv(RUN/'NONTARGET_SELECTION_COUNTS.csv',index=False)
    save_npz(RUN/'data/BALANCED_MAPS.npz',cr_safe=cr_safe,cr_absolute=cr_abs,all_safe=all_safe,all_absolute=all_abs,all_response_censored_safe=all_censored_safe,participant_ids=PIDS,channels=CHANNELS,frequencies_hz=FREQ,times_ms=(POST-600)/FS*1000,codes=CODES)
    (RUN/'logs/balanced_maps.json').write_text(json.dumps({'status':'PASS','participant_n':16,'common_code_count':10,'specifications':['frequency_safe_ersp','absolute_log_post'],'deterministic':'equal unweighted code means before feature discovery'},indent=2)+'\n')
if __name__=='__main__':main()
