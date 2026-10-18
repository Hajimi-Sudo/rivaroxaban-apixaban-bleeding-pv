from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT = Path('/root/autodl-tmp/Frontiers_in_Pharmacology/upgrade_2026-10-14/04_manuscript/latex/figures')
RES = Path('/root/autodl-tmp/Frontiers_in_Pharmacology/upgrade_2026-10-14/results')
OUT.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8.5, 'axes.titlesize':9,
                     'axes.labelsize':8.5, 'xtick.labelsize':8, 'ytick.labelsize':8,
                     'figure.dpi':160, 'savefig.dpi':300, 'axes.linewidth':0.7,
                     'pdf.fonttype':42, 'ps.fonttype':42})
BLUE = '#0072B2'
ORANGE = '#D55E00'
GREY = '#4D4D4D'
LIGHT_GREY = '#F2F2F2'

def save(fig, name):
    fig.savefig(OUT / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUT / f'{name}.png', bbox_inches='tight', dpi=300)
    plt.close(fig)

# Figure 1: audited endpoint-specific FAERS forest plot.
csv = RES / 'streaming_2015Q1_2026Q1_final' / 'outcome_group_counts.csv'
df = pd.read_csv(csv)
labels = {
    'menstrual_specific':'Menstrual-specific',
    'abnormal_uterine':'Abnormal uterine bleeding',
    'uterine_haemorrhage':'Uterine haemorrhage',
    'vaginal_haemorrhage':'Vaginal haemorrhage',
    'historical_broad_composite':'Historical broad composite',
}
df['label'] = df['outcome_group'].map(labels)
df = df.iloc[[0,1,2,3,4]].copy()
y = np.arange(len(df))[::-1]
fig, ax = plt.subplots(figsize=(7.2,4.4))
for yi, row in zip(y, df.itertuples()):
    color = BLUE if row.comparator_event_stable_20 else ORANGE
    ax.plot([row.ror_ci95_lower,row.ror_ci95_upper],[yi,yi],color=color,lw=1.8,solid_capstyle='round')
    ax.scatter(row.ror,[yi],s=34,color=color,zorder=3,edgecolor='white',linewidth=0.5)
    if not row.comparator_event_stable_20:
        ax.text(min(row.ror_ci95_upper*1.15, 55), yi, 'sparse', va='center', fontsize=7.5, color=ORANGE)
ax.axvline(1,color=GREY,lw=0.8,ls=(0,(3,2)))
ax.set_xscale('log'); ax.set_yticks(y); ax.set_yticklabels(df['label'])
ax.set_xlabel('Reporting odds ratio (log scale; 95% CI)')
ax.set_title('')
ax.grid(axis='x',color='#D9D9D9',lw=0.45); ax.spines[['top','right']].set_visible(False)
ax.set_xlim(0.2,60)
ax.tick_params(length=3,width=0.6)
save(fig,'figure1_faers_endpoint_forest')

# Figure 2: transparent flow from reports to mutually exclusive cohorts.
fig, ax = plt.subplots(figsize=(8.0,4.8)); ax.axis('off'); ax.set_xlim(0,10); ax.set_ylim(0,6.2)
boxes = [
    (0.2,2.4,1.9,1.2,'FAERS quarterly files\n2015 Q1–2026 Q1\n15,508,934\ndeduplicated cases'),
    (2.55,2.4,1.8,1.2,'Female reports\n7,990,400'),
    (4.8,2.4,1.9,1.2,'Female age 15–49\n7,299 eligible reports'),
    (7.25,3.8,2.2,1.0,'Rivaroxaban-only\n4,196 reports'),
    (7.25,2.35,2.2,1.0,'Apixaban-only\n1,986 reports'),
    (7.25,0.9,2.2,1.0,'Outside active\ncontrast\n1,117 reports'),
]
for x,y0,w,h,t in boxes:
    ax.add_patch(FancyBboxPatch((x,y0),w,h,boxstyle='round,pad=0.04,rounding_size=0.04',fc=LIGHT_GREY,ec=GREY,lw=0.8))
    ax.text(x+w/2,y0+h/2,t,ha='center',va='center',fontsize=7.8,color='#222222')
arrows=[((2.15,3.0),(2.48,3.0)),((4.4,3.0),(4.73,3.0)),((6.75,3.15),(7.18,4.25)),((6.75,3.0),(7.18,2.85)),((6.75,2.85),(7.18,1.4))]
for a,b in arrows: ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=9,lw=0.85,color=GREY))
ax.text(5,0.25,'Dual-exposure reports excluded: 0; remaining reports were not used in the active-comparator contrast.',ha='center',fontsize=7.5,color=GREY)
save(fig,'figure2_faers_flow')

# Figure 3: 45-quarter audited baseline vs 46-quarter temporal update.
# A distinct three-panel result figure uses an annotated heatmap, a change
# intensity tile map, and a donut chart; no lines or bars are used here.
base = {'Menstrual-specific':(0.5228857841,0.3452612691,0.7918917286),
        'Uterine haemorrhage':(7.9863129184,4.2077543229,15.1580128344),
        'Vaginal haemorrhage':(6.3433728448,4.5737030140,8.7977682252),
        'Historical broad composite':(3.9218569592,3.1283188546,4.9166861576)}
upd = {'Menstrual-specific':(0.523,0.352,0.778),
       'Uterine haemorrhage':(8.244,4.344,15.647),
       'Vaginal haemorrhage':(6.548,4.722,9.080),
       'Historical broad composite':(3.858,3.093,4.811)}
labs=list(base.keys())
fig=plt.figure(figsize=(8.2,6.2))
gs=fig.add_gridspec(2,2,height_ratios=[1.35,1.0],hspace=0.56,wspace=0.42)
ax=fig.add_subplot(gs[0,:]); ax2=fig.add_subplot(gs[1,0]); ax3=fig.add_subplot(gs[1,1])
# Panel A: annotated heatmap of the audited and updated ROR estimates.
matrix=np.array([[base[lab][0],upd[lab][0]] for lab in labs])
im=ax.imshow(np.log10(matrix),cmap='YlGnBu',aspect='auto',vmin=-0.35,vmax=1.05)
ax.set_xticks([0,1]); ax.set_xticklabels(['45-quarter audited','46-quarter update'])
ax.set_yticks(np.arange(len(labs))); ax.set_yticklabels(['Menstrual-specific','Uterine haemorrhage','Vaginal haemorrhage','Broad composite'])
ax.set_title('(A) Endpoint-specific ROR matrix',loc='left',fontsize=9,pad=7)
ax.set_xlabel('Temporal window'); ax.set_ylabel('Endpoint')
for i in range(matrix.shape[0]):
    for j in range(matrix.shape[1]):
        ax.text(j,i,f'{matrix[i,j]:.2f}',ha='center',va='center',fontsize=8,
                color='white' if np.log10(matrix[i,j])>0.45 else '#222222',fontweight='bold')
cb=fig.colorbar(im,ax=ax,fraction=0.025,pad=0.02)
cb.set_label('log10(ROR)',fontsize=7.5); cb.ax.tick_params(labelsize=7)
ax.set_xticks(np.arange(-.5,2,1),minor=True); ax.set_yticks(np.arange(-.5,len(labs),1),minor=True)
ax.grid(which='minor',color='white',lw=1.2); ax.tick_params(which='minor',bottom=False,left=False)
ax.spines[['top','right']].set_visible(False)

# Panel B: change-intensity tiles with exact percentage annotations.
delta=np.array([upd[lab][0]/base[lab][0]-1 for lab in labs])*100
lim=max(abs(delta).max(),4.0)
im2=ax2.imshow(np.array([delta]),cmap='RdBu_r',aspect='auto',vmin=-lim,vmax=lim)
ax2.set_xticks(np.arange(len(labs))); ax2.set_xticklabels(['Menstrual','Uterine','Vaginal','Broad'],rotation=25,ha='right')
ax2.set_yticks([0]); ax2.set_yticklabels(['ROR change'])
ax2.set_title('(B) Point-estimate change after 2026 Q2',loc='left',fontsize=9,pad=7)
for j,val in enumerate(delta):
    ax2.text(j,0,f'{val:+.1f}%',ha='center',va='center',fontsize=9,fontweight='bold',color='#222222')
ax2.set_xticks(np.arange(-.5,len(labs),1),minor=True); ax2.set_yticks([-0.5,0.5],minor=True)
ax2.grid(which='minor',color='white',lw=1.2); ax2.tick_params(which='minor',bottom=False,left=False)
ax2.set_xlabel('Endpoint'); ax2.spines[['top','right']].set_visible(False)
cb2=fig.colorbar(im2,ax=ax2,fraction=0.05,pad=0.06)
cb2.set_label('Change (%)',fontsize=7.5); cb2.ax.tick_params(labelsize=7)

# Panel C: donut chart of the 119 newly added reports.
addition_labels=['Rivaroxaban-only','Apixaban-only','Outside contrast']
addition_values=[30,79,10]
addition_colors=[BLUE,ORANGE,GREY]
donut_center=(0.28,0.50)
wedges,_=ax3.pie(addition_values,colors=addition_colors,startangle=90,counterclock=False,
                 center=donut_center,radius=0.30,
                 wedgeprops=dict(width=0.11,edgecolor='white',linewidth=1.2))
ax3.text(donut_center[0],donut_center[1]+0.02,'119',ha='center',va='center',fontsize=13,fontweight='bold',color='#222222')
ax3.text(donut_center[0],donut_center[1]-0.10,'new reports',ha='center',va='center',fontsize=6.5,color=GREY)
ax3.set_xlim(0,1); ax3.set_ylim(0,1); ax3.set_aspect('equal'); ax3.axis('off')
ax3.set_title('(C) Composition of the temporal increment',loc='left',fontsize=9,pad=7)
ax3.legend(wedges,[f'{lab}  (n={val}; {val/sum(addition_values)*100:.1f}%)' for lab,val in zip(addition_labels,addition_values)],
           frameon=False,fontsize=6.8,loc='center left',bbox_to_anchor=(0.55,0.50),ncol=1,
           handlelength=1.0,handletextpad=0.45)
for panel in (ax,ax2,ax3): panel.tick_params(length=3,width=0.6)
fig.subplots_adjust(top=0.91,bottom=0.23,left=0.14,right=0.98)
save(fig,'figure3_temporal_update')

# Figure 4: contextual spontaneous-report estimates across databases.
# A three-panel result figure shows the database-specific RORs, the report
# event fractions, and the raw event-cell sizes that determine sparse-cell
# interpretation. These are not a harmonized or pooled endpoint.
context_rows=[
 ('FAERS\n45-quarter',3.9218569592,3.1283188546,4.9166861576,BLUE),
 ('Canada\nVigilance',4.45,1.63,15.22,ORANGE),
 ('JADER',1.85,0.78,4.42,GREY),
]
y=np.arange(len(context_rows))[::-1]
fig=plt.figure(figsize=(8.2,6.2))
gs=fig.add_gridspec(2,2,height_ratios=[1.35,1.0],hspace=0.50,wspace=0.38)
ax=fig.add_subplot(gs[0,:]); ax2=fig.add_subplot(gs[1,0]); ax3=fig.add_subplot(gs[1,1])
# Panel A: confidence-capsule plot. Each translucent capsule spans the 95% CI
# on a log ROR axis, while the diamond marks the point estimate.
xx=np.arange(len(context_rows))
ax.axhspan(0.9,1.1,color='#E8E8E8',alpha=0.65,zorder=0)
for xi,(lab,est,lo,hi,col) in zip(xx,context_rows):
    ax.add_patch(plt.Rectangle((xi-0.18,lo),0.36,hi-lo,facecolor=col,edgecolor=col,
                               alpha=0.22,lw=1.2,zorder=1))
    ax.scatter(xi,est,s=82,marker='D',color=col,edgecolor='white',linewidth=1.0,zorder=3)
    ax.text(xi,hi*1.18,f'{est:.2f}',ha='center',va='bottom',fontsize=8.0,color=col,fontweight='bold')
    ax.text(xi,lo/1.16,f'{lo:.2f}–{hi:.2f}',ha='center',va='top',fontsize=6.8,color=GREY)
ax.set_yscale('log'); ax.set_ylim(0.55,25)
ax.set_xticks(xx); ax.set_xticklabels(['FAERS\n45-quarter','Canada\nVigilance','JADER'])
ax.set_ylabel('Database-specific contextual ROR (log scale)')
ax.set_title('(A) Contextual ROR capsules (95% CI)',loc='left',fontsize=9,pad=7)
ax.set_yticks([1,2,5,10,20]); ax.set_yticklabels(['1','2','5','10','20'])
ax.grid(axis='y',color='#D9D9D9',lw=0.45); ax.spines[['top','right']].set_visible(False)

db_labels=['FAERS','Canada','JADER']
riva_events=np.array([665,47,18]); riva_non=np.array([3531,103,42])
apix_events=np.array([91,5,15]); apix_non=np.array([1895,49,65])
rr=riva_events/(riva_events+riva_non)*100
aa=apix_events/(apix_events+apix_non)*100
xx=np.arange(len(db_labels)); bw=0.32
ax2.bar(xx-bw/2,rr,width=bw,color=BLUE,label='Rivaroxaban-only')
ax2.bar(xx+bw/2,aa,width=bw,color=ORANGE,label='Apixaban-only')
ax2.set_xticks(xx); ax2.set_xticklabels(db_labels)
ax2.set_ylabel('Event reports (%)')
ax2.set_title('(B) Event fraction within each database-specific cohort',loc='left',fontsize=8.6,pad=7)
ax2.grid(axis='y',color='#D9D9D9',lw=0.45); ax2.spines[['top','right']].set_visible(False)
ax2.legend(frameon=False,fontsize=7.0,loc='upper center',bbox_to_anchor=(0.5,-0.20),
           ncol=2,columnspacing=0.7,handletextpad=0.4)
for xi,v in zip(xx-bw/2,rr): ax2.text(xi,v+1.2,f'{v:.1f}',ha='center',fontsize=7.2,color=BLUE)
for xi,v in zip(xx+bw/2,aa): ax2.text(xi,v+1.2,f'{v:.1f}',ha='center',fontsize=7.2,color=ORANGE)

ax3.scatter(np.zeros(3),y,s=np.sqrt(riva_events)*18,color=BLUE,alpha=0.9,label='Rivaroxaban event cell')
ax3.scatter(np.ones(3),y,s=np.sqrt(apix_events)*18,color=ORANGE,alpha=0.9,label='Apixaban event cell')
for yi,rv,av in zip(y,riva_events,apix_events):
    ax3.text(-0.10,yi,str(rv),ha='right',va='center',fontsize=7.5,color=BLUE)
    ax3.text(1.10,yi,str(av),ha='left',va='center',fontsize=7.5,color=ORANGE)
ax3.set_xlim(-0.45,1.45); ax3.set_xticks([0,1]); ax3.set_xticklabels(['Rivaroxaban','Apixaban'],rotation=18,ha='right')
ax3.set_yticks(y); ax3.set_yticklabels(db_labels)
ax3.set_ylabel('Database'); ax3.set_title('(C) Raw event-cell size (n)',loc='left',fontsize=9,pad=7)
ax3.grid(axis='x',color='#D9D9D9',lw=0.45); ax3.spines[['top','right']].set_visible(False)
ax3.legend(frameon=False,fontsize=6.8,loc='upper center',bbox_to_anchor=(0.5,-0.30),
           ncol=2,columnspacing=0.7,handletextpad=0.35,markerscale=0.55)
for panel in (ax,ax2,ax3): panel.tick_params(length=3,width=0.6)
fig.subplots_adjust(top=0.91,bottom=0.36,left=0.14,right=0.98)
save(fig,'figure5_contextual_database_forest')

# Figure 5: endpoint-specific event burden by active comparator.
# This complements the ROR forest plot with directly interpretable event rates
# and uses a grouped horizontal bar layout rather than another interval plot.
rate_df = df.copy()
rate_df['rivaroxaban_rate_per_1000'] = rate_df['a_rivaroxaban_event'] / (rate_df['a_rivaroxaban_event'] + rate_df['b_rivaroxaban_non_event']) * 1000
rate_df['apixaban_rate_per_1000'] = rate_df['c_apixaban_event'] / (rate_df['c_apixaban_event'] + rate_df['d_apixaban_non_event']) * 1000
rate_df = rate_df.iloc[[0,1,2,3,4]].copy()
order = np.arange(len(rate_df))[::-1]
fig, ax = plt.subplots(figsize=(7.2,4.6))
bar_h = 0.32
ax.barh(order + bar_h/2, rate_df['rivaroxaban_rate_per_1000'], height=bar_h,
        color=BLUE, label='Rivaroxaban-only')
ax.barh(order - bar_h/2, rate_df['apixaban_rate_per_1000'], height=bar_h,
        color=ORANGE, label='Apixaban-only')
ax.set_yticks(order); ax.set_yticklabels(rate_df['label'])
ax.set_xlabel('Endpoint reports per 1,000 active-comparator reports')
ax.grid(axis='x',color='#D9D9D9',lw=0.45); ax.spines[['top','right']].set_visible(False)
ax.legend(frameon=False,loc='upper left',bbox_to_anchor=(0.56,1.0),fontsize=7.5,ncol=2)
for yi, rv, av in zip(order, rate_df['rivaroxaban_rate_per_1000'], rate_df['apixaban_rate_per_1000']):
    xmax = max(rv, av)
    ax.text(xmax+0.35, yi+bar_h/2, f'{rv:.1f}', va='center', fontsize=7.2, color=BLUE)
    ax.text(xmax+0.35, yi-bar_h/2, f'{av:.1f}', va='center', fontsize=7.2, color=ORANGE)
save(fig,'figure5_endpoint_event_burden')

print('Generated:', ', '.join(p.name for p in OUT.glob('figure*.pdf')))
