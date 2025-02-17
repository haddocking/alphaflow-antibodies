import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from pathlib import Path

pdb_files_low = ['7bnv', '7dk2', '7f7e', '7kez', '7kql',
                 '7msq', '7n4j', '7np1', '7nx3', '7phu',
                 '7phw', '7pr0', '7ps2', '7ps6', '7qny',
                 '7rfb', '7s11']

rmsd_data_flow = pd.read_csv("../data/alphaflow_1000models_rmsd.csv", sep="\t")
rmsd_data_fold = pd.read_csv("../data/colabfold_25models_rmsd.csv", sep="\t")

# violinplot with the best values
lowest_h3_rmsd_fold = []
lowest_h3_rmsd_flow = []

lowest_50_h3_rmsd_flow = []
lowest_100_h3_rmsd_flow = []
lowest_200_h3_rmsd_flow = []
lowest_500_h3_rmsd_flow = []

for pdb in pdb_files_low:
    pdb_rmsd_data_fold = rmsd_data_fold[rmsd_data_fold["pdb"] == pdb]
    pdb_rmsd_data_flow = rmsd_data_flow[rmsd_data_flow["pdb"] == pdb]
    # lowest h3 rmsd
    lowest_h3_fold = pdb_rmsd_data_fold.loc[pdb_rmsd_data_fold["rmsd_h3"] == pdb_rmsd_data_fold["rmsd_h3"].min()]
    lowest_h3_flow = pdb_rmsd_data_flow.loc[pdb_rmsd_data_flow["rmsd_h3"] == pdb_rmsd_data_flow["rmsd_h3"].min()]
    print(f"pdb: {pdb}, lowest h3-rmsd af2: {lowest_h3_fold['rmsd_h3'].values[0]} ({lowest_h3_fold['model'].values[0]}), lowest h3-rmsd afl: {lowest_h3_flow['rmsd_h3'].values[0]} ({lowest_h3_flow['model'].values[0]})")
    lowest_h3_rmsd_fold.append(lowest_h3_fold["rmsd_h3"].values[0])
    lowest_h3_rmsd_flow.append(lowest_h3_flow["rmsd_h3"].values[0])
    # first 50 models (model-id < 50)
    first_50_h3_flow = pdb_rmsd_data_flow[pdb_rmsd_data_flow["model-id"] < 50]
    lowest_50_h3_rmsd_flow.append(first_50_h3_flow["rmsd_h3"].min())
    # first 100 models (model-id < 100)
    first_100_h3_flow = pdb_rmsd_data_flow[pdb_rmsd_data_flow["model-id"] < 100]
    lowest_100_h3_rmsd_flow.append(first_100_h3_flow["rmsd_h3"].min())
    # first 200 models (model-id < 200)
    first_200_h3_flow = pdb_rmsd_data_flow[pdb_rmsd_data_flow["model-id"] < 200]
    lowest_200_h3_rmsd_flow.append(first_200_h3_flow["rmsd_h3"].min())
    # first 500 models (model-id < 500)
    first_500_h3_flow = pdb_rmsd_data_flow[pdb_rmsd_data_flow["model-id"] < 500]
    lowest_500_h3_rmsd_flow.append(first_500_h3_flow["rmsd_h3"].min())

print(f"lowest_500_h3_rmsd_flow = {lowest_500_h3_rmsd_flow}")

# violinplot with all the N data
plt.violinplot([lowest_h3_rmsd_fold, lowest_h3_rmsd_flow, lowest_500_h3_rmsd_flow, 
                lowest_200_h3_rmsd_flow, lowest_100_h3_rmsd_flow, lowest_50_h3_rmsd_flow])
plt.xticks([1, 2, 3, 4, 5, 6], ["AF2", "AFL\n(N=1000)", "AFL\n(N=500)", 
                                "AFL\n(N=200)", "AFL\n(N=100)", "AFL\n(N=50)"])
plt.ylabel("H3 RMSD [$\AA$]")
plt.title("Lowest H3 RMSD")
plt.savefig(Path("figures", "alphaflow", "lowest_h3_rmsd_violinplot_N.png"), dpi=300)
plt.close()

#
# now let's put the two plots in one, with the scatter plot in an inset.
#
fig, ax = plt.subplots(figsize=(5, 6))
p = ax.violinplot([lowest_h3_rmsd_fold, lowest_h3_rmsd_flow], showextrema=False)
p["bodies"][0].set_facecolor('blue')
p["bodies"][0].set_alpha(1)
p["bodies"][1].set_facecolor('orange')
p["bodies"][1].set_alpha(0.9)

ax.set_xticks([1, 2])
ax.set_xticklabels(["AF2", "AFL"], fontsize=14)
ax.set_ylabel("H3 RMSD [$\AA$]", fontsize=14)
ax.set_title("Lowest H3 RMSD", fontsize=16)
ax_inset = fig.add_axes([0.5, 0.533, 0.47, 0.4])
ax_inset.scatter(lowest_h3_rmsd_fold, lowest_h3_rmsd_flow, s=5, color="black")
ax_inset.set_xlabel("AF2 H3 RMSD [$\AA$]")
ax_inset.set_ylabel("AFL H3 RMSD [$\AA$]")
# place the label closer to the axis
ax_inset.yaxis.set_label_coords(-0.125, 0.5)
ax_inset.xaxis.set_label_coords(0.5, -0.125)

max_rmsd = max(max(lowest_h3_rmsd_fold), max(lowest_h3_rmsd_flow))
max_y_int = round(max_rmsd + 0.5) 
max_y = max_y_int + 0.25
ax_inset.set_xlim(0, max_y)
#ax_inset.set_ylim(0, 3.5)
ax_inset.set_ylim(0, max_y)
ticks = [n for n in range(0, max_y_int+1)]
ax_inset.set_xticks(ticks)
ax_inset.set_yticks(ticks)
# change fontsize of the ticks
ax_inset.tick_params(axis='both', which='major', labelsize=8)

# add vertical line at 3.0
ax_inset.axvline(x=3.0, color="gray", linestyle="--", linewidth=0.75)
# add diagonal
ax_inset.plot([0, max_y], [0, max_y], color="gray", linestyle="solid", linewidth=0.75, alpha=0.5)
# add LOW-80-W text on the right and LOW-80-R on the left
ax_inset.text(4.45, 0.1, "LOW-80-W", fontsize=10)
ax_inset.text(0.2, 0.1, "LOW-80-R", fontsize=10)
plt.tight_layout()
plt.savefig(Path("figures", "alphaflow", "lowest_h3_rmsd_violinplot_inset.png"), dpi=400)
plt.close()

# now do the histogram of the improvements
improvements = []
impr_high3 = [] # improvements when lowest_h3_rmsd_fold > 3.0
impr_low3 = [] # improvements when lowest_h3_rmsd_fold < 3.0
pdb_files_bad = []
for i in range(len(pdb_files_low)):
    improvements.append(lowest_h3_rmsd_fold[i] - lowest_h3_rmsd_flow[i])
    if lowest_h3_rmsd_fold[i] > 3.0:
        print(f"pdb: {pdb_files_low[i]}")
        pdb_files_bad.append(pdb_files_low[i])
        impr_high3.append(lowest_h3_rmsd_fold[i] - lowest_h3_rmsd_flow[i])
    else:
        impr_low3.append(lowest_h3_rmsd_fold[i] - lowest_h3_rmsd_flow[i])

print(f"avg improvement: {np.mean(improvements):.2f} +- {np.std(improvements):.2f}")
print(f"avg improvement for high h3 rmsd: {np.mean(impr_high3):.2f} +- {np.std(impr_high3):.2f} ({len(impr_high3)} entries)")
print(f"avg improvement for low h3 rmsd: {np.mean(impr_low3):.2f} +- {np.std(impr_low3):.2f} ({len(impr_low3)} entries)")

# now the avg improvements for the first 50, 100, 200, 500 models
print(f"avg improvement for the first 50 models: {np.mean(lowest_h3_rmsd_fold) - np.mean(lowest_50_h3_rmsd_flow):.2f}")
print(f"avg improvement for the first 100 models: {np.mean(lowest_h3_rmsd_fold) - np.mean(lowest_100_h3_rmsd_flow):.2f}")
print(f"avg improvement for the first 200 models: {np.mean(lowest_h3_rmsd_fold) - np.mean(lowest_200_h3_rmsd_flow):.2f}")
print(f"avg improvement for the first 500 models: {np.mean(lowest_h3_rmsd_fold) - np.mean(lowest_500_h3_rmsd_flow):.2f}")

print("\nclustering analysis\n")
print(f"pdb_files_bad (LOW-80-W dataset) = {pdb_files_bad}")

improvement_df = pd.read_csv("../data/improvement_data.csv", sep="\t")
cl_methods = improvement_df["cl_method"].unique()
# drop - from cl_methods
cl_methods = cl_methods[cl_methods != "-"]

# print the unique pdbs
unique_pdbs = improvement_df["pdb"].unique()
assert (unique_pdbs == pdb_files_low).all(), "The unique pdbs are not the same as the expected ones"

# one plot per cl_method
for cl_method in cl_methods:
    unique_models = improvement_df[improvement_df["cl_method"] == cl_method]["cl_num"].unique()
    for cl_num in unique_models:
        sub_df = improvement_df[(improvement_df["cl_method"] == cl_method) & (improvement_df["cl_num"] == cl_num)]
        avg_improvement = sub_df.groupby("topn")["impr"].mean().reset_index()
        print(f"cl_num {cl_num} ({cl_method}): avg_improvement {[(avg_improvement.iloc[i, 0], avg_improvement.round(2).iloc[i, 1]) for i in range(len(avg_improvement))]}")
        plt.scatter(avg_improvement["topn"], avg_improvement["impr"], label=cl_num)
        plt.plot(avg_improvement["topn"], avg_improvement["impr"], linewidth = 0.2, linestyle="-")
    # now the unclustered data
    sub_df = improvement_df[(improvement_df["cl_method"] == "-") & (improvement_df["cl_num"] == "ALL")]
    avg_improvement = sub_df.groupby("topn")["impr"].mean().reset_index()
    plt.scatter(avg_improvement["topn"], avg_improvement["impr"], label="ALL")
    plt.plot(avg_improvement["topn"], avg_improvement["impr"], linewidth = 0.1)
    # nice legend
    plt.xscale("log")
    plt.ylabel("Avg improvement vs best AF2 model\n(H3 rmsd) $[\AA]$", fontsize=16)
    plt.xlabel("Top N Alphaflow models considered", fontsize=16)
    # would like to put proper labels on the x axis
    plt.xticks([10, 20, 50, 100, 200, 500, 1000], ["10", "20", "50", "100", "200", "500", "1000"], fontsize=16)
    plt.legend(title="Number of clusters generated", loc="upper left", ncol=2)
    plt.tight_layout()
    plt.savefig(Path("figures", "alphaflow", f"avg_improvement_vs_topn_{cl_method}.png"), dpi=300)
    plt.close()

print(f"\n\n\nNow the bad pdb files ({pdb_files_bad})\n\n\n")
# now we consider only the bad pdbs
improvement_df_bad = improvement_df[improvement_df["pdb"].isin(pdb_files_bad)]

# one plot per cl_method
for cl_method in cl_methods:
    unique_models = improvement_df_bad[improvement_df_bad["cl_method"] == cl_method]["cl_num"].unique()
    for cl_num in unique_models:
        sub_df = improvement_df_bad[(improvement_df_bad["cl_method"] == cl_method) & (improvement_df_bad["cl_num"] == cl_num)]
        avg_improvement = sub_df.groupby("topn")["impr"].mean().reset_index()
        # print(f"cl_num {cl_num} ({cl_method}): avg_improvement {avg_improvement.round(2).to_dict()}")
        print(f"cl_num {cl_num} ({cl_method}): avg_improvement {[(avg_improvement.iloc[i, 0], avg_improvement.round(2).iloc[i, 1]) for i in range(len(avg_improvement))]}")
        plt.scatter(avg_improvement["topn"], avg_improvement["impr"], label=cl_num)
        plt.plot(avg_improvement["topn"], avg_improvement["impr"], linewidth = 0.4, linestyle="--")
    # now the unclustered data
    sub_df = improvement_df_bad[(improvement_df_bad["cl_method"] == "-") & (improvement_df_bad["cl_num"] == "ALL")]
    avg_improvement = sub_df.groupby("topn")["impr"].mean().reset_index()
    plt.scatter(avg_improvement["topn"], avg_improvement["impr"], label="ALL")
    plt.plot(avg_improvement["topn"], avg_improvement["impr"], linewidth = 0.4, linestyle="--")
    # nice legend
    plt.xscale("log")
    plt.ylabel("Avg improvement vs best AF2 model\n(H3 rmsd) $[\AA]$", fontsize=16)
    plt.xlabel("Top N Alphaflow models considered", fontsize=16)
    # would like to put proper labels on the x axis
    plt.xticks([10, 20, 50, 100, 200, 500, 1000], ["10", "20", "50", "100", "200", "500", "1000"], fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(title="Number of clusters generated", loc="upper left", ncol=2)
    plt.tight_layout()
    plt.savefig(Path("figures", "alphaflow", f"avg_improvement_vs_topn_{cl_method}_bad.png"), dpi=300)
    plt.close()