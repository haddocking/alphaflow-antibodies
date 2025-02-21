import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr

def make_loop_violinplot(data, title, output_name):
    print(f"making violinplot {output_name}")
    fig, ax = plt.subplots()
    plots = ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
    for i, pc in enumerate(plots['bodies']):
        if i < 3:
            pc.set_facecolor(f"orange")
        else:
            pc.set_facecolor(f"skyblue")
    #ax.set_title(title)
    ax.set_ylabel("RMSD [$\AA$]", fontsize=20)
    ax.set_xticks([1, 2, 3, 4, 5, 6])
    ax.set_xticklabels(["H1", "H2", "H3", "L1", "L2", "L3"], fontsize=20)
    ## Set the color of the violin patches
    averages = [np.mean(el) for el in data]
    for i, mean in enumerate(averages):
        ax.text(i + 1, mean, f'${mean:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=15)
    ## change the fontsize of the y labels
    ax.tick_params(axis='y', labelsize=15)
    plt.savefig(Path("figures", "alphafold", output_name), dpi=300)
    plt.close()

rmsd_data = pd.read_csv("../data/colabfold_25models_rmsd.csv", sep="\t")
plddt_data = pd.read_csv("../data/colabfold_25models_plddt.csv", sep="\t")

# first plot: best ranked model: loop accuracy
best_ranked_loop_rmsd = rmsd_data.loc[rmsd_data["model"].str.contains("rank_001")]
# how many entries have a value of rmsd_h3 > 3.0
print(f"number of entries with rmsd_h3 > 3.0: {best_ranked_loop_rmsd[best_ranked_loop_rmsd['rmsd_h3'] > 3.0].shape[0]}")
# print those entries
print(f"entries with rmsd_h3 > 3.0: {best_ranked_loop_rmsd[best_ranked_loop_rmsd['rmsd_h3'] > 3.0]['pdb']} (in % {best_ranked_loop_rmsd[best_ranked_loop_rmsd['rmsd_h3'] > 3.0].shape[0] / best_ranked_loop_rmsd.shape[0] * 100:.2f})")
data = [[], [], [], [], [], []]
for i, row in best_ranked_loop_rmsd.iterrows():
    data[0].append(row["rmsd_h1"])
    data[1].append(row["rmsd_h2"])
    data[2].append(row["rmsd_h3"])
    data[3].append(row["rmsd_l1"])
    data[4].append(row["rmsd_l2"])
    data[5].append(row["rmsd_l3"])
make_loop_violinplot(data, "Best ranked alphafold model", "rank001_loop_rmsd_violinplot.png")

# SI FIG.1 : for each loop we select the RMSD corresponding to the highest pLDDT
highest_plddt_data = [[], [], [], [], [], []]
for n, pdb in enumerate(plddt_data["pdb"].unique()):
    for m, loop in enumerate(["h1", "h2", "h3", "l1", "l2", "l3"]):
        pdb_rmsd_data = rmsd_data[rmsd_data["pdb"] == pdb]
        pdb_plddt_data = plddt_data[plddt_data["pdb"] == pdb]
        best_plddt = pdb_plddt_data.loc[pdb_plddt_data[f"plddt_{loop}"] == pdb_plddt_data[f"plddt_{loop}"].max()]
        best_rmsd = pdb_rmsd_data.loc[pdb_rmsd_data["model"] == best_plddt["model"].values[0]]
        highest_plddt_data[m].append(best_rmsd[f"rmsd_{loop}"].values[0])
        # print(f"pdb: {pdb}, loop: {loop}, plddt {best_plddt[f'plddt_{loop}']} rmsd: {best_rmsd[f'rmsd_{loop}'].values[0]}")
make_loop_violinplot(highest_plddt_data, "Highest pLDDT loops", "highest_plddt_loop_rmsd_violinplot.png")

# second plot: best ranked h3 loop (according to h3 plddt) vs h3 rmsd
best_ranked_h3plddt = []
best_ranked_h3rmsd = []
best_h3plddt = []
best_h3rmsd = []

pdb_files_low = []

for pdb in rmsd_data["pdb"].unique():
    pdb_rmsd_data = rmsd_data[rmsd_data["pdb"] == pdb]
    pdb_plddt_data = plddt_data[plddt_data["pdb"] == pdb]
    # highest plddt for h3
    best_ranked_h3 = pdb_plddt_data.loc[pdb_plddt_data["plddt_h3"] == pdb_plddt_data["plddt_h3"].max()]
    # if the highest plddt is < 80 append this pdb to pdb_files_low
    if best_ranked_h3["plddt_h3"].values[0] < 80:
        pdb_files_low.append(pdb)
    # corresponding rmsd
    best_ranked_h3_rmsd = pdb_rmsd_data.loc[pdb_rmsd_data["model"] == best_ranked_h3["model"].values[0]]
    # best antibody according to the alphafold ranking scheme
    # this means select the model that contains the string "rank_001"
    best_ranked_antibody = pdb_plddt_data[pdb_plddt_data["model"].str.contains("rank_001")]
    # now the corresponding rmsd
    best_ranked_antibody_rmsd = pdb_rmsd_data.loc[pdb_rmsd_data["model"] == best_ranked_antibody["model"].values[0]]
    
    best_ranked_h3plddt.append(best_ranked_h3["plddt_h3"].values[0])
    best_ranked_h3rmsd.append(best_ranked_h3_rmsd["rmsd_h3"].values[0])
    best_h3plddt.append(best_ranked_antibody["plddt_h3"].values[0])
    best_h3rmsd.append(best_ranked_antibody_rmsd["rmsd_h3"].values[0])

print(f"pdb files with low plddt: {pdb_files_low} (len = {len(pdb_files_low)})")
# calculate how many of these pdbs have a best rmsd > 3.0
pdb_files_low_high_rmsd = []
for pdb in pdb_files_low:
    pdb_rmsd_data = rmsd_data[rmsd_data["pdb"] == pdb]
    pdb_rmsd_data.sort_values("rmsd_h3", inplace=True)
    print(f"pdb: {pdb} rmsd_h3: {pdb_rmsd_data['rmsd_h3'].values[0]} corr. {pdb_rmsd_data['model'].values[0]}")
    if pdb_rmsd_data["rmsd_h3"].values[0] > 3.0:
        pdb_files_low_high_rmsd.append(pdb)
print(f"pdb files with low plddt and high rmsd: {pdb_files_low_high_rmsd} (len = {len(pdb_files_low_high_rmsd)})")

plt.scatter(best_ranked_h3plddt, best_ranked_h3rmsd)
plt.ylabel("RMSD [$\AA$]")
plt.xlabel("H3 plddt")
plt.title("H3 RMSD vs H3 plddt")
# add the correlation somewhere on the plot
correlation = pearsonr(best_ranked_h3plddt, best_ranked_h3rmsd)[0]
print(f"best h3: h3 rmsd-plddt correlation : {correlation:.2f}")
#plt.text(90.0, 5.0, f"correlation: {pearsonr(best_ranked_h3rmsd, best_ranked_h3plddt)[0]:.2f}", fontsize=12, ha='center')
plt.savefig(Path("figures", "alphafold", "H3rmsd_vs_H3plddt_besth3plddt.png"), dpi=300)
plt.close()

# now the same plot but coloured by length
df_loop = pd.read_csv("../data/loop_data.csv", sep="\t")
h3_lengths = []
for pdb in rmsd_data["pdb"].unique():
    h3_loop = df_loop[df_loop["PDB"] == pdb.lower()]["H3"].values[0]
    h3_lengths.append(len(h3_loop))

#colormap = plt.cm.get_cmap('PuOr')
#colors = colormap(scaled_z)
plt.scatter(best_ranked_h3plddt, best_ranked_h3rmsd, c=h3_lengths, cmap="cividis")
#plt.scatter(best_ranked_h3plddt, best_ranked_h3rmsd, c=h3_lengths)
cbar = plt.colorbar()
# add legend to colorbar
cbar.set_label('Number of residues', loc="center", rotation=90, fontsize=12)
plt.ylabel("RMSD [$\AA$]")
plt.xlabel("H3 plddt")
plt.title("H3 RMSD vs H3 plddt")
plt.savefig(Path("figures", "alphafold", "H3rmsd_vs_H3plddt_besth3plddt_coloured.png"), dpi=300)
plt.close()

# 30/1/2025: we create a unique plot for figure 1 using 

fig, ax = plt.subplots(1,2, figsize=(12,6), width_ratios=[3, 4])
plots = ax[0].violinplot(data, showmeans=False, showmedians=False, showextrema=False)
for i, pc in enumerate(plots['bodies']):
    if i < 3:
        pc.set_facecolor(f"orange")
    else:
        pc.set_facecolor(f"skyblue")
#ax.set_title(title)
ax[0].set_ylabel("RMSD [$\AA$]", fontsize=18)
ax[0].set_xticks([1, 2, 3, 4, 5, 6])
ax[0].set_xticklabels(["H1", "H2", "H3", "L1", "L2", "L3"], fontsize=16)
ax[0].set_xlabel("Loop", fontsize=18)
## Set the color of the violin patches
averages = [np.mean(el) for el in data]
for i, mean in enumerate(averages):
    ax[0].text(i + 1, mean, f'${mean:.2f}\ \AA$', ha='center', va='bottom', color='black', fontsize=16)
## change the fontsize of the y labels
ax[0].tick_params(axis='y', labelsize=16)
ax[1].scatter(best_ranked_h3plddt, best_ranked_h3rmsd, c=h3_lengths, cmap="cividis")
# ax[1].set_ylabel("RMSD [$\AA$]")
ax[1].set_xlabel("H3 pLDDT", fontsize=18)
cbar = plt.colorbar(ax[1].scatter(best_ranked_h3plddt, best_ranked_h3rmsd, c=h3_lengths, cmap="cividis"), ax=ax[1])
cbar.set_label('Number of residues', loc="center", rotation=90, fontsize=16)
# cbar ticks fontsize
cbar.ax.tick_params(labelsize=16)
# x labels : starting from 65 to 100 spaced by 5
ax[1].set_xticks(np.arange(65, 101, 5))
ax[1].tick_params(axis='x', labelsize=16)
# no y ticks
ax[1].set_yticks([])
# adding a and b labels 
# ax[0].text(-1.5, 105, "a)", fontsize = 18)
# similar to this but using relative numbers
ax[0].text(-0.15, 0.93, "a)", fontsize = 24, transform=ax[0].transAxes) 
ax[1].text(0.01, 0.93, "b)", fontsize = 24, transform=ax[1].transAxes)
plt.subplots_adjust(hspace=0, wspace=0)
plt.tight_layout()
plt.savefig(Path("figures", "alphafold", "figure1.png"), dpi=400)
plt.close()

# now with the best ranked antibody
plt.scatter(best_h3plddt, best_h3rmsd)
plt.ylabel("RMSD [$\AA$]")
plt.xlabel("H3 plddt")
plt.title("H3 RMSD vs H3 plddt")
# add the correlation somewhere on the plot
correlation = pearsonr(best_h3rmsd, best_h3plddt)[0]
print(f"best antibody: h3 rmsd-plddt correlation : {correlation:.2f}")
#plt.text(90.0, 5.0, f"correlation: {pearsonr(best_h3rmsd, best_h3plddt)[0]:.2f}", fontsize=12, ha='center')
plt.savefig(Path("figures", "H3rmsd_vs_H3plddt_bestantibody.png"), dpi=300)
plt.close()

# analysis between runs with or without antigen
antigen_rmsd_data = pd.read_csv("../data/colabfold_25models_antigen_rmsd.csv", sep="\t")
antigen_plddt_data = pd.read_csv("../data/colabfold_25models_antigen_plddt.csv", sep="\t")
unique_pdbs = antigen_rmsd_data["pdb"].unique()
print(f"\nComparison between AF2 with and without antigen ({len(unique_pdbs)} pdbs) \n")
print(f"unique_pdbs: {unique_pdbs}")
non_in_low = [pdb for pdb in unique_pdbs if pdb not in pdb_files_low]
print(f"pdbs in antigen data but not in pdb_files_low: {non_in_low}")
not_in_antigen = [pdb for pdb in pdb_files_low if pdb not in unique_pdbs]
print(f"pdbs in pdb_files_low but not in antigen data: {not_in_antigen}")
assert len(unique_pdbs) == len(pdb_files_low), f"antigen data has {len(unique_pdbs)} unique pdbs, while pdb_files_low has {len(pdb_files_low)} unique pdbs"
# assessment: must be done only for the pdb_files_low
comparison_data = []
for pdb in pdb_files_low:
    pdb_rmsd_data = rmsd_data[rmsd_data["pdb"] == pdb]
    pdb_plddt_data = plddt_data[plddt_data["pdb"] == pdb]
    pdb_antigen_rmsd_data = antigen_rmsd_data[antigen_rmsd_data["pdb"] == pdb]
    pdb_antigen_plddt_data = antigen_plddt_data[antigen_plddt_data["pdb"] == pdb]
    # highest plddt for h3
    best_ranked_h3 = pdb_plddt_data.loc[pdb_plddt_data["plddt_h3"] == pdb_plddt_data["plddt_h3"].max()]
    # corresponding rmsd
    best_ranked_h3_rmsd = pdb_rmsd_data.loc[pdb_rmsd_data["model"] == best_ranked_h3["model"].values[0]]
    # best antibody according to the alphafold ranking scheme
    # this means select the model that contains the string "rank_001"
    best_ranked_antibody = pdb_plddt_data[pdb_plddt_data["model"].str.contains("rank_001")]
    # now the corresponding rmsd
    best_ranked_antibody_rmsd = pdb_rmsd_data.loc[pdb_rmsd_data["model"] == best_ranked_antibody["model"].values[0]]
    # antigen

    best_ranked_antigen_h3 = pdb_antigen_plddt_data.loc[pdb_antigen_plddt_data["plddt_h3"] == pdb_antigen_plddt_data["plddt_h3"].max()]
    # print(f"best_ranked_antigen_h3: {best_ranked_antigen_h3}")
    best_ranked_antigen_h3_rmsd = pdb_antigen_rmsd_data.loc[pdb_antigen_rmsd_data["model"] == best_ranked_antigen_h3["model"].values[0]]
    best_ranked_antigen = pdb_antigen_plddt_data[pdb_antigen_plddt_data["model"].str.contains("rank_001")]
    best_ranked_antigen_rmsd = pdb_antigen_rmsd_data.loc[pdb_antigen_rmsd_data["model"] == best_ranked_antigen["model"].values[0]]
    # overall best model
    antigen_best = pdb_antigen_rmsd_data.loc[pdb_antigen_rmsd_data["rmsd_h3"] == pdb_antigen_rmsd_data["rmsd_h3"].min()]
    antibody_best = pdb_rmsd_data.loc[pdb_rmsd_data["rmsd_h3"] == pdb_rmsd_data["rmsd_h3"].min()]
    # now the comparison
    comparison_data.append([pdb, best_ranked_h3_rmsd['rmsd_h3'].values[0], best_ranked_antigen_h3_rmsd['rmsd_h3'].values[0],
                            best_ranked_antibody_rmsd['rmsd_h3'].values[0], best_ranked_antigen_rmsd['rmsd_h3'].values[0],
                            antibody_best['rmsd_h3'].values[0], antigen_best['rmsd_h3'].values[0],
                            best_ranked_h3['plddt_h3'].values[0], best_ranked_antigen_h3['plddt_h3'].values[0]])

# do some stats on comparison_data: what is the average difference between the best model with and without antigen
# and the best ranked model with and without antigen
# and the best ranked h3 with and without antigen
print("\nStats on comparison data\n")
print(f"overall best model antibody only has an average H3-rmsd of {np.mean([el[5] for el in comparison_data]):.2f} +- {np.std([el[5] for el in comparison_data]):.2f}")
print(f"overall best model antigen only has an average H3-rmsd of {np.mean([el[6] for el in comparison_data]):.2f} +- {np.std([el[6] for el in comparison_data]):.2f}")
print(f"best ranked model antibody only has an average H3-rmsd of {np.mean([el[3] for el in comparison_data]):.2f}")
print(f"best ranked model antigen only has an average H3-rmsd of {np.mean([el[4] for el in comparison_data]):.2f}")
print(f"best ranked h3 antibody only has an average H3-rmsd of {np.mean([el[1] for el in comparison_data]):.2f}")
print(f"best ranked h3 antigen only has an average H3-rmsd of {np.mean([el[2] for el in comparison_data]):.2f}")
print(f"best ranked h3 antibody only has an average H3-plddt of {np.mean([el[7] for el in comparison_data]):.2f}")
print(f"best ranked h3 antigen only has an average H3-plddt of {np.mean([el[8] for el in comparison_data]):.2f}")

comparison_df = pd.DataFrame(comparison_data, columns=["pdb", "best_ranked_h3", "best_ranked_antigen_h3", "best_ranked_antibody", "best_ranked_antigen", "antibody_best", "antigen_best", "best_ranked_h3_plddt", "best_ranked_antigen_plddt"])
# create a boxplot for the comparison data
fig, ax = plt.subplots()
plots = ax.boxplot([comparison_df["best_ranked_h3"], comparison_df["best_ranked_antigen_h3"], comparison_df["best_ranked_antibody"], comparison_df["best_ranked_antigen"], comparison_df["antibody_best"], comparison_df["antigen_best"]], showmeans=False, showfliers=False)
ax.set_title("Comparison of H3 RMSD values")
ax.set_ylabel("RMSD [$\AA$]")
ax.set_xticks([1, 2, 3, 4, 5, 6])
ax.set_xticklabels(["BR H3-AB", "BR H3-ANT", "BR AB", "BR ANT", "AB best", "ANT best"])
# put the xlabels vertically
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(Path("figures", "alphafold", "comparison_antigen.png"), dpi=300)