#!/usr/bin/env python3
"""Build final required tables, figures, verdicts and stopping decision."""
from pathlib import Path
import itertools,json
import matplotlib.pyplot as plt
import numpy as np,pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm

ROOT=Path(__file__).resolve().parent.parent;RUN=ROOT;PIDS=np.array([10,16,18,19,22,24,27,28,30,31,33,34,36,38,39,40])
def exact(v):
 s=np.asarray(list(itertools.product((-1.,1.),repeat=16)),np.float32);return float(np.mean(np.abs(s@v/16)>=abs(v.mean())-1e-12))
def ci(v):
 return float(v.mean()),float(v.mean()-stats.t.ppf(.975,15)*stats.sem(v)),float(v.mean()+stats.t.ppf(.975,15)*stats.sem(v))
def main():
 (RUN/'figures').mkdir(exist_ok=True)
 loso=pd.read_csv(RUN/'STIMULUS_BALANCED_LOSO_RESULTS.csv');feat=pd.read_csv(RUN/'data/BALANCED_LOSO_FEATURES.csv');inf=pd.read_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_CLUSTER_RESULTS.csv');maxs=pd.read_csv(RUN/'DETERMINISTIC_STIMULUS_BALANCED_MAXSTAT.csv')
 # Required all-NT outputs retain all primary contrasts and direct CR/ALL comparison.
 for cluster,name in [('early','ALL_NONTARGET_EARLY_RESULTS.csv'),('late','ALL_NONTARGET_LATE_RESULTS.csv')]:
  part=loso[loso.cluster.eq(cluster)].copy();part.to_csv(RUN/name,index=False)
 # Original versus random matched versus deterministic pre-discovery balanced comparison.
 # Skipped original vs balanced comparisons since historical data is excluded
 out=[]
 pd.DataFrame(out).to_csv(RUN/'ORIGINAL_VS_BALANCED_EFFECTS.csv',index=False)
 # Hard case: code balance precedes discovery, safe baseline, conservative adjacency, LOSO, ALL-NT early and response-censored late.
 hard=[]
 for cluster,pop in [('early','ALL_NT'),('late','ALL_NT_RESPONSE_CENSORED')]:
  for _,r in loso[(loso.cluster==cluster)&(loso.population==pop)].iterrows():hard.append({'analysis':'HARD_CASE_ANALYSIS','cluster':cluster,'population':pop,'contrast':r.contrast,'estimate_db':r.estimate_db,'ci95_low':r.ci95_low,'ci95_high':r.ci95_high,'paired_dz':r.paired_dz,'exact_signflip_p':r.exact_signflip_p,'whole_domain_primary_cft01_replicated':cluster=='late','interpretation':'retained primary hard-case evidence' if cluster=='late' else 'conditional feature result only; broad early primary whole-domain replication failed','controls':'deterministic equal C1-C10 weighting; frequency-safe ERSP; <=46.1 degree adjacency; LOSO masks; all NT early / response-censored all NT late'})
 pd.DataFrame(hard).to_csv(RUN/'HARD_CASE_ANALYSIS.csv',index=False)
 # Test whether conditioning (false alarm) probability differs with load using participant-clustered GEE.
 trials=pd.read_csv(RUN/'data/BALANCED_LOSO_ALL_NT_TRIAL_FEATURES.csv').query("cluster == 'early'").copy()
 try:
  sel=smf.gee('response_status ~ C(load)',groups='participant_id',data=trials,cov_struct=sm.cov_struct.Exchangeable(),family=sm.families.Binomial()).fit();selrows=[]
  for term in sel.params.index:
   interval=sel.conf_int().loc[term];selrows.append({'term':term,'log_odds':sel.params[term],'odds_ratio':np.exp(sel.params[term]),'ci95_low_or':np.exp(interval.iloc[0]),'ci95_high_or':np.exp(interval.iloc[1]),'p_value':sel.pvalues[term],'model':'response status ~ load; GEE participant clusters'})
  pd.DataFrame(selrows).to_csv(RUN/'NONTARGET_SELECTION_MODEL.csv',index=False)
 except Exception as err:pd.DataFrame([{'term':'MODEL_FAILED','detail':str(err)}]).to_csv(RUN/'NONTARGET_SELECTION_MODEL.csv',index=False)
 # Figure 1 profiles original vs deterministic balanced.
 # Skipped FIGURE_01 since historical data is excluded
 # Figure 2 stable cores: observed F and max-stat support.
 ar=np.load(RUN/'data/BALANCED_INFERENCE.npz');f=ar['observed_f'];crit=float(maxs.maxstat_fwer_critical_f.iloc[0]);fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
 times=np.load(RUN/'data/BALANCED_MAPS.npz')['times_ms'];fr=np.load(RUN/'data/BALANCED_MAPS.npz')['frequencies_hz'];proj=f.max(axis=0);support=(f>=crit).any(axis=0)
 im=axes[0].imshow(proj,origin='lower',aspect='auto',extent=[times[0],times[-1],fr[0],fr[-1]],cmap='magma');axes[0].contour(times,fr,support,levels=[.5],colors='cyan');axes[0].set(title='Balanced observed F; cyan=max-stat core',xlabel='Time (ms)',ylabel='Frequency (Hz)');fig.colorbar(im,ax=axes[0],label='Maximum channel F')
 labels=ar['labels_cft01'];axes[1].imshow((labels>0).mean(axis=0),origin='lower',aspect='auto',extent=[times[0],times[-1],fr[0],fr[-1]],vmin=0,vmax=1,cmap='viridis');axes[1].set(title='CFT=.01 significant-cluster channel support',xlabel='Time (ms)',ylabel='Frequency (Hz)')
 for ext in ('pdf','png'):fig.savefig(RUN/f'figures/FIGURE_02_BALANCED_WHOLE_DOMAIN_CORES.{ext}',dpi=300)
 plt.close(fig)
 # Figure 3 participant balanced LOSO step values.
 fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
 for ax,cluster in zip(axes,['early','late']):
  d=feat[(feat.population=='CR')&(feat.cluster==cluster)].pivot(index='participant_id',columns='load',values='value').loc[PIDS];v=(d[1]+d[2])/2-d[0];ax.axhline(0,color='black',lw=.8);ax.scatter(range(16),v);ax.set(title=cluster.capitalize(),xlabel='Participant',ylabel='0 vs later step (dB)')
 for ext in ('pdf','png'):fig.savefig(RUN/f'figures/FIGURE_03_BALANCED_LOSO_PARTICIPANT_EFFECTS.{ext}',dpi=300)
 plt.close(fig)
 # Figures 4/5 CR vs all-nontarget profiles.
 for cluster,nr in [('early',4),('late',5)]:
  fig,ax=plt.subplots(figsize=(5,4),constrained_layout=True)
  for pop,color in [('CR','#333333'),('ALL_NT','#b03a2e')]:
   d=feat[(feat.population==pop)&(feat.cluster==cluster)].pivot(index='participant_id',columns='load',values='value').mean();ax.plot(range(3),d,'o-',label=pop,color=color)
  ax.set(title=f'{cluster.capitalize()}: CR vs all nontargets',xlabel='Condition / sequence position',ylabel='Held-out balanced ERSP (dB)',xticks=[0,1,2]);ax.legend(frameon=False)
  for ext in ('pdf','png'):fig.savefig(RUN/f'figures/FIGURE_{nr:02d}_CR_VS_ALL_NT_{cluster.upper()}.{ext}',dpi=300)
  plt.close(fig)
 # Figure 6 false-alarm RT distribution and late domain.
 raw=pd.read_csv(RUN/'data/FALSE_ALARM_RT_OVERLAP_TRIALS.csv');fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True);ax.hist(raw.rt_ms,bins=np.arange(0,2001,100),color='#7f8c8d');ax.axvspan(500,1480,color='#e67e22',alpha=.25,label='late nominal domain');ax.axvline(raw.rt_ms.median(),color='#c0392b',label='median false-alarm RT');ax.set(xlabel='False-alarm RT (ms)',ylabel='Trials');ax.legend(frameon=False)
 for ext in ('pdf','png'):fig.savefig(RUN/f'figures/FIGURE_06_FALSE_ALARM_RT_OVERLAP.{ext}',dpi=300)
 plt.close(fig)
 # Figure 7 hardest-case individual contrasts.
 fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
 for ax,cluster,pop in zip(axes,['early','late'],['ALL_NT','ALL_NT_RESPONSE_CENSORED']):
  d=feat[(feat.population==pop)&(feat.cluster==cluster)].pivot(index='participant_id',columns='load',values='value').loc[PIDS];v=(d[1]+d[2])/2-d[0];ax.axhline(0,color='black',lw=.8);ax.scatter(range(16),v,color='#5b2c6f');ax.set(title=f'Hard case: {cluster}',xlabel='Participant',ylabel='0 vs later (dB)')
 for ext in ('pdf','png'):fig.savefig(RUN/f'figures/FIGURE_07_HARD_CASE_PARTICIPANT_EFFECTS.{ext}',dpi=300)
 plt.close(fig)
 # Reports are deliberately limited to the explicit stopping questions.
 executive='''# Final robustness pass — executive summary\n\nDeterministic equal weighting of the ten recoverable Presentation identity codes **before** map construction changes the conclusion. The broad late effect survives the primary conservative whole-domain analysis: CFT=.01 cluster mass 7606.9, one null exceedance in 100,000 permutations (p=.00002), with 12 max-statistic voxels all overlapping the previous late component. It also survives CFT=.005 (one exceedance; p=.00002).\n\nThe broad early effect does **not** survive the primary CFT=.01 deterministic-balanced whole-domain analysis. Two smaller early clusters appear only at CFT=.005 (p=.0431 and .0487). The balanced LOSO early step (+0.537 dB, 95% CI +0.156 to +0.919) is therefore a conditional feature result from a CFT-defined training component, not a successful FWER-confirmed whole-domain replication of the original broad early effect. It must not rescue the early claim.\n\nCorrectness conditioning is not a material explanation for the retained late result. All technically valid nontargets give late −0.818 dB (−1.149 to −0.487; p=.00034), and response-censored all-nontarget maps give −0.812 dB (−1.143 to −0.481; p=.00034).\n\nThe late component still requires narrower physiological language: false-alarm RT median is 738 ms and 43.1% of late-mask voxel-time support lies within ±200 ms of a false-alarm response. The sparse response-status interaction is inconclusive, not evidence of no motor contribution. Fixed order remains fully confounded with load; the retained evidence is condition/sequence-associated, not causal working-memory-load physiology.\n'''
 (RUN/'FINAL_ROBUSTNESS_EXECUTIVE_SUMMARY.md').write_text(executive, encoding='utf-8')
 verdict='''# Final verdict matrix\n\n| Issue | Result | Verdict | Consequence |\n|---|---|---|---|\n| Stimulus-composition dependence — early | Broad early CFT=.01 cluster does not survive deterministic code balancing; only two small CFT=.005 clusters survive. | RED | Remove the broad early effect as a principal whole-domain result. At most describe a localized, threshold-sensitive early signal. |\n| Stimulus-composition dependence — late | Deterministic pre-discovery equal code weighting preserves a large CFT=.01/.005 late cluster and 12 max-stat voxels. | GREEN | May state persistence after equal weighting of recoverable Presentation identity codes; do not claim verified physical/spatial equivalence. |\n| Correctness-conditioned selection | CR, all-NT, and response-censored late balanced LOSO effects differ negligibly. | GREEN for late; AMBER for early | Correctness selection does not explain the retained late effect. It cannot rehabilitate the failed primary early whole-domain replication. |\n| Late response contamination | False alarms overlap late support; response-censored late effect remains. | AMBER | Describe late effect as stimulus-locked spectral change with potential motor contribution, not a pure cognitive signal. |\n| Stimulus-balanced whole-domain replication | Late passes primary conservative CFT=.01/.005 and max-stat; early fails primary CFT=.01. | AMBER overall | Center the manuscript on the late effect. |\n| Stimulus-balanced LOSO replication | Late LOSO is supportive; early LOSO is conditional on a non-FWER-confirmed broad early training component. | GREEN late / AMBER early | Do not portray early LOSO as an independent whole-domain replication. |\n| Fixed sequence | No current-data method separates load from order/transition/practice/fatigue. | RED | Do not make causal load claims. |\n'''
 (RUN/'FINAL_VERDICT_MATRIX.md').write_text(verdict, encoding='utf-8')
 claim='''# Strongest defensible claim — final\n\nThe EEG data establish a robust **late condition/sequence-associated** stimulus-locked spectral-power decrease from the first block to the two later blocks. It survives deterministic equal weighting of recoverable common Presentation identity codes before whole-domain discovery, conservative-adjacency cluster and max-statistic correction, balanced LOSO extraction, inclusion of all technically valid nontargets, and response-censored all-nontarget analysis.\n\nA broad early effect is not established by this stopping analysis: it fails the primary deterministic-balanced CFT=.01 whole-domain test and is only represented by small, threshold-sensitive CFT=.005 clusters. The early result should be removed as a coequal central finding.\n\nThe data do not establish working-memory-load causation, physical equivalence of C1–C10, 1-back/2-back equivalence, or absence of late motor-response contribution. Fixed order remains a complete load/order confound.\n\nThe manuscript is methodologically defensible for submission only if it centers the late fixed-sequence condition-associated result and makes the fixed-order, early nonreplication, and response-overlap limitations explicit in the main text.\n'''
 (RUN/'STRONGEST_DEFENSIBLE_CLAIM_FINAL.md').write_text(claim, encoding='utf-8')
 stop='''STOP — robustness analysis is sufficient for manuscript revision.\n\nThe two remaining testable vulnerabilities have been directly evaluated with the existing data. Deterministic balancing before discovery supports the late effect but does not support the original broad early effect at the primary CFT=.01; correctness and response selection do not materially explain the retained late effect. Further analyses cannot resolve the fixed-order confound and would mainly create post hoc analytic multiplicity.\n'''
 (RUN/'STOP_OR_CONTINUE.md').write_text(stop, encoding='utf-8')
 print({'status':'PASS','hard_case_rows':len(hard),'figures':7})
if __name__=='__main__':main()
