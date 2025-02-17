import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

def extract_bars(df, n_list, npdbs):
    """
    Extract the success rate for each n in n_list
    """
    bar_acc = []
    bar_med = []
    bar_high = []
    bar_near_acc = []
    for n in n_list:
        sr_acc = df.loc[(df["n"] == n) & (df["tx_acc"] != 0)].shape[0]/npdbs
        bar_acc.append(sr_acc)
        sr_med = df.loc[(df["n"] == n) & (df["tx_med"] != 0)].shape[0]/npdbs
        bar_med.append(sr_med)
        sr_high = df.loc[(df["n"] == n) & (df["tx_high"] != 0)].shape[0]/npdbs
        bar_high.append(sr_high)        
    return bar_acc, bar_med, bar_high

ns = [1,5,10,20,50,100]
tns = [f"T{n}" for n in ns]
ns_aranged = np.arange(len(ns))

df_pe_emref_b = pd.read_csv(Path("../data", "docking_sr", "pe_emref_b.tsv"),sep="\t")
df_pe_emref_u = pd.read_csv(Path("../data", "docking_sr", "pe_emref_u.tsv"),sep="\t")
df_pe_afens_emref_b = pd.read_csv(Path("../data", "docking_sr", "afens_pe_emref_b.tsv"),sep="\t")
df_pe_afens_emref_u = pd.read_csv(Path("../data", "docking_sr", "afens_pe_emref_u.tsv"),sep="\t")
df_epivag_emref_b = pd.read_csv(Path("../data", "docking_sr", "emref_b.tsv"),sep="\t")
df_epivag_emref_u = pd.read_csv(Path("../data", "docking_sr", "emref_u.tsv"),sep="\t")
df_epivag_afens_emref_b = pd.read_csv(Path("../data", "docking_sr", "afens_emref_b.tsv"),sep="\t")
df_epivag_afens_emref_u = pd.read_csv(Path("../data", "docking_sr", "afens_emref_u.tsv"),sep="\t")

NPDBS = df_pe_emref_b["pdb"].nunique()

bars_pe_emref_b_acc, bars_pe_emref_b_med, bars_pe_emref_b_high = extract_bars(df_pe_emref_b, ns, NPDBS)
bars_pe_emref_u_acc, bars_pe_emref_u_med, bars_pe_emref_u_high = extract_bars(df_pe_emref_u, ns, NPDBS)
bars_afens_pe_emref_b_acc, bars_afens_pe_emref_b_med, bars_afens_pe_emref_b_high = extract_bars(df_pe_afens_emref_b, ns, NPDBS)
bars_afens_pe_emref_u_acc, bars_afens_pe_emref_u_med, bars_afens_pe_emref_u_high = extract_bars(df_pe_afens_emref_u, ns, NPDBS)
# epivag
bars_epivag_emref_b_acc, bars_epivag_emref_b_med, bars_epivag_emref_b_high = extract_bars(df_epivag_emref_b, ns, NPDBS)
bars_epivag_emref_u_acc, bars_epivag_emref_u_med, bars_epivag_emref_u_high = extract_bars(df_epivag_emref_u, ns, NPDBS)
bars_afens_epivag_emref_b_acc, bars_afens_epivag_emref_b_med, bars_afens_epivag_emref_b_high = extract_bars(df_epivag_afens_emref_b, ns, NPDBS)
bars_afens_epivag_emref_u_acc, bars_afens_epivag_emref_u_med, bars_afens_epivag_emref_u_high = extract_bars(df_epivag_afens_emref_u, ns, NPDBS)

# now let's make a whole figure (2x2)
fig, ax = plt.subplots(2,2, figsize=(12,12))
width = 0.35

# pe
ax[0,0].bar(ns_aranged-width/2, bars_pe_emref_b_acc, width, label='Acceptable AFL', color='#a5cee2')
ax[0,0].bar(ns_aranged-width/2, bars_pe_emref_b_med, width, label='Medium AFL', color='#b2df8a')
ax[0,0].bar(ns_aranged-width/2, bars_pe_emref_b_high, width, label='High AFL', color='#33a02b')
# now let's compare to afens emref_b data (bars_afens_emref_b_acc) slightly shifted to the right
ax[0,0].bar(ns_aranged+width/2, bars_afens_pe_emref_b_acc, width, label='Acceptable AF2', color='bisque')
ax[0,0].bar(ns_aranged+width/2, bars_afens_pe_emref_b_med, width, label='Medium AF2', color='sandybrown')
# ax[0,0].bar(ns_aranged+width/2, bars_afens_pe_emref_b_high, width, label='High AF2', color='darkred')

ax[0,1].bar(ns_aranged-width/2, bars_pe_emref_u_acc, width, color='#a5cee2')
ax[0,1].bar(ns_aranged-width/2, bars_pe_emref_u_med, width, color='#b2df8a')
ax[0,1].bar(ns_aranged-width/2, bars_pe_emref_u_high, width, color='#33a02b')
# now let's compare to afens emref_u data (bars_afens_emref_u_acc) slightly shifted to the right
ax[0,1].bar(ns_aranged+width/2, bars_afens_pe_emref_u_acc, width, color='bisque')
ax[0,1].bar(ns_aranged+width/2, bars_afens_pe_emref_u_med, width, color='sandybrown')
ax[0,1].bar(ns_aranged+width/2, bars_afens_pe_emref_u_high, width, color='darkorange')

ax[1,0].bar(ns_aranged-width/2, bars_epivag_emref_b_acc, width, color='#a5cee2')
ax[1,0].bar(ns_aranged-width/2, bars_epivag_emref_b_med, width, color='#b2df8a')
ax[1,0].bar(ns_aranged-width/2, bars_epivag_emref_b_high, width, color='#33a02b')

ax[1,0].bar(ns_aranged+width/2, bars_afens_epivag_emref_b_acc, width, color='bisque')
ax[1,0].bar(ns_aranged+width/2, bars_afens_epivag_emref_b_med, width, color='sandybrown')
ax[1,0].bar(ns_aranged+width/2, bars_afens_epivag_emref_b_high, width, color='darkorange')

ax[1,1].bar(ns_aranged-width/2, bars_epivag_emref_u_acc, width, color='#a5cee2')
ax[1,1].bar(ns_aranged-width/2, bars_epivag_emref_u_med, width, color='#b2df8a')
ax[1,1].bar(ns_aranged-width/2, bars_epivag_emref_u_high, width, color='#33a02b')

ax[1,1].bar(ns_aranged+width/2, bars_afens_epivag_emref_u_acc, width, color='bisque')
ax[1,1].bar(ns_aranged+width/2, bars_afens_epivag_emref_u_med, width, color='sandybrown')
ax[1,1].bar(ns_aranged+width/2, bars_afens_epivag_emref_u_high, width, color='darkorange')
        
for i in range(2):
    ax[0,i].set_ylim(0,1.01)
    ax[1,i].set_ylim(0,1.01)
    ax[i,0].set_ylabel('Success Rate (%)', fontsize=18)
    ax[1,i].set_xticks(ns_aranged)
    ax[1,i].set_xticklabels(tns, fontsize=14)
    # for ax[0,i] delete xticks
    ax[0,i].set_xticks([])
    ax[i,0].set_yticks(np.arange(0,1.01,0.2))
    ax[i,0].set_yticklabels([f"{el}" for el in np.arange(0,101,20)], fontsize=14)
    # for ax[i,1] delete yticks
    ax[i,1].set_yticks([])
    # horizontal bars at 20-40-60-80
    for j in range(1,5):
        ax[i,0].axhline(j*0.2, color='grey', linestyle='dashed', linewidth=1)
        ax[i,1].axhline(j*0.2, color='grey', linestyle='dashed', linewidth=1)

ax[0,0].set_title('Bound antigen', fontsize=18)
ax[0,1].set_title('Unbound antigen', fontsize=18)
# put a y label on the right side of ax[0,1]
ax[0,1].set_ylabel('Para-Epi', fontsize=18)
ax[0,1].yaxis.set_label_position("right")
ax[1,1].set_ylabel('CDR-EpiVag', fontsize=18)
ax[1,1].yaxis.set_label_position("right")
# I would like these two labels to be a bit more to the right
ax[0,1].yaxis.set_label_coords(1.03,0.5)
ax[1,1].yaxis.set_label_coords(1.03,0.5)

fig.legend(loc='lower center', bbox_to_anchor=(0.5, -0.055), ncol=3, fontsize=14, bbox_transform=fig.transFigure)

plt.tight_layout()
plt.savefig("figures/docking_table.png", bbox_inches='tight', dpi=400)