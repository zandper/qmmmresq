#!/usr/bin/env python3
"""
Create DataFrame with residues as rows, frames as columns of lambda_max values.
"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from schrodinger import structure

# Residue classification
pos_charged = {"ARG", "LYS", "HIS"}
neg_charged = {"ASP", "GLU", "ASH"}
polar = {"SER", "THR", "ASN", "GLN", "TYR", "CYS"}
nonpolar = {"ALA", "VAL", "LEU", "ILE", "MET", "PRO", "PHE", "TRP", "GLY"}

def classify_residue(rescode):
    rescode = str(rescode).upper()
    if rescode in pos_charged:
        return "(+) Charge"
    elif rescode in neg_charged:
        return "(-) Charge"
    elif rescode in polar:
        return "Polar"
    elif rescode in nonpolar:
        return "Nonpolar"
    else:
        return "Other/Nonprotein"

# Color map (matching your style)
COLOR_MAP = {
    "(+) Charge": "#4a90e2",      # muted blue
    "(-) Charge": "#d64541",      # soft red
    "Polar": "#9b59b6",           # muted purple
    "Nonpolar": "#f1c40f",        # softer yellow
    "Other/Nonprotein": "#95a5a6" # grayish soft
}

def create_frame_contrib_matrix(param_folders, mae_st):
    """Create matrix of contributions with residue codes."""
    
    # Build mapping from (molnum, resnum) to rescode
    rescode_map = {}
    for residue in mae_st.residue:
        key = (residue.molecule_number, residue.resnum)
        rescode_map[key] = residue.pdbres.strip()
    
    rows = []
    for folder in param_folders:
        parts = os.path.basename(folder).split('_')
        frame_num = parts[-2]
        
        for txt_file in glob.glob(os.path.join(folder, "*.txt")):
            with open(txt_file) as f:
                for line in f:
                    parts_line = line.strip().split()
                    if len(parts_line) >= 5:
                        molnum = int(parts_line[0])
                        resnum = int(parts_line[1])
                        rows.append({
                            'molnum': molnum,
                            'resnum': resnum,
                            'rescode': rescode_map.get((molnum, resnum), 'UNK'),
                            'frame': frame_num,
                            'contrib': float(parts_line[-2])
                        })
    
    if not rows:
        return pd.DataFrame()
    
    df = pd.DataFrame(rows)
    
    # Pivot to get frames as columns
    matrix = df.pivot(index=['molnum', 'resnum', 'rescode'], columns='frame', values='contrib')
    matrix.columns.name = 'frame'
    
    # Calculate mean and std across frames
    matrix['mean'] = matrix.mean(axis=1)
    matrix['std'] = matrix.std(axis=1)
    matrix['n_frames'] = matrix.count(axis=1)
    
    return matrix


def plot_bar_grouped(df, fname, title="", top_n=10, bar_height=0.8):
    """
    Horizontal bar plot colored by residue class.
    Highlights top N residues by absolute mean contribution in bold on y-axis.
    """
    # Prepare data
    plot_df = df.reset_index().copy()
    
    # Add classification and labels
    plot_df['res_class'] = plot_df['rescode'].apply(classify_residue)
    plot_df['res_label'] = plot_df['rescode'] + plot_df['resnum'].astype(str)
    
    # Sort by class and mean contribution
    plot_df = plot_df.sort_values(by=["res_class", "mean"], ascending=[True, True]).reset_index(drop=True)
    
    # Identify top N residues by |mean|
    top_indices = plot_df["mean"].abs().nlargest(top_n).index
    
    # Build plot
    fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.4)))
    
    y_positions = np.arange(len(plot_df))
    colors = plot_df["res_class"].map(COLOR_MAP)
    
    # Add error bars if std available
    ax.barh(
        y_positions,
        plot_df["mean"],
        xerr=plot_df.get('std', None),
        color=colors,
        edgecolor="black",
        height=bar_height,
        capsize=3,
        error_kw={'elinewidth': 1.5, 'capthick': 1.5}
    )
    
    # Y labels: bold top N
    y_labels = []
    for i, label in enumerate(plot_df["res_label"]):
        if i in top_indices:
            y_labels.append(f"$\\bf{{{label}}}$")
        else:
            y_labels.append(label)
    
    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=14)
    ax.set_xlabel(r'$\Delta \lambda$ Contribution (nm)', fontsize=16)
    ax.set_ylabel("", fontsize=0)
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    # Add zero line
    ax.axvline(x=0, color='black', linestyle='-', linewidth=1, alpha=0.5)
    
    ax.invert_yaxis()
    
    # Legend
    used_classes = plot_df["res_class"].unique()
    handles = [
        plt.Line2D([0], [0], color=COLOR_MAP[k], lw=8)
        for k in used_classes
    ]
    legend = ax.legend(handles, used_classes, title="", fontsize=12, loc="lower right")
    legend.get_frame().set_linewidth(1.5)
    
    # Tick params
    ax.tick_params(axis='x', labelsize=12, width=1.5, length=6)
    ax.tick_params(axis='y', labelsize=12, width=1.5, length=6)
    
    plt.tight_layout()
    outname = f"{fname}_grouped_bar.png"
    print(f"Saving matplotlib grouped bar as {outname}")
    plt.savefig(outname, dpi=300, bbox_inches='tight')
    plt.close()
    
    return plot_df


if __name__ == "__main__":
    mae_path = './000001_geopt.01.mae'
    mae_st = structure.StructureReader.read(mae_path)
    
    param_folders = sorted(glob.glob('./parameters_*'))
    matrix = create_frame_contrib_matrix(param_folders, mae_st)
    
    # Plot all residues
    plot_bar_grouped(matrix, 'all_residues', title='All Residue Contributions', top_n=10)
    
    # Plot top 25 by absolute contribution
    top25_indices = matrix['mean'].abs().nlargest(25).index
    matrix_top25 = matrix.loc[top25_indices]
    plot_bar_grouped(matrix_top25, 'top25_residues', title='Top 25 Residue Contributions', top_n=5)
    
    # Save data
    matrix.reset_index().to_csv('contributions.csv', index=False)
    print(f"Processed {len(matrix)} residues")
    print("Saved: contributions.csv")