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
import sys
from scipy import stats
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
        frame_num = parts[1]
        
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
                            'contrib': float(parts_line[-2])#*-1 #REMOVE THIS LATER
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
    matrix['min'] = matrix.min(axis=1)
    matrix['max'] = matrix.max(axis=1)

    # Calculate confidence interval
    ci_level=0.95
    # CI = mean ± t * (std / sqrt(n))
    z_score = stats.t.ppf((1 + ci_level) / 2, df=matrix['n_frames'] - 1)
    matrix['ci_half_width'] = z_score * (matrix['std'] / np.sqrt(matrix['n_frames']))
    matrix['ci_lower'] = matrix['mean'] - matrix['ci_half_width']
    matrix['ci_upper'] = matrix['mean'] + matrix['ci_half_width']
    return matrix


def plot_bar_grouped(df, fname, title="", top_n=10, bar_height=0.8, x_min=None, x_max=None):
    plot_df = df.reset_index().copy()
    plot_df['res_class'] = plot_df['rescode'].apply(classify_residue)
    plot_df['res_label'] = plot_df['rescode'] + plot_df['resnum'].astype(str)
    
    # Define custom order: (+) Charge, (-) Charge, Polar, Nonpolar, Other/Nonprotein
    class_order = ["(+) Charge", "(-) Charge", "Polar", "Nonpolar", "Other/Nonprotein"]
    
    # Convert res_class to categorical with specified order
    plot_df['res_class'] = pd.Categorical(plot_df['res_class'], categories=class_order, ordered=True)
    
    # Sort by res_class (categorical order) then by mean
    plot_df = plot_df.sort_values(by=["res_class", "mean"], ascending=[True, True]).reset_index(drop=True)
    
    top_indices = plot_df["mean"].abs().nlargest(top_n).index

    # Identify frame columns (everything that isn't a summary stat or metadata)
    meta_cols = {'molnum', 'resnum', 'rescode', 'res_class', 'res_label',
                 'mean', 'std', 'n_frames', 'min', 'max'}
    frame_cols = [c for c in plot_df.columns if c not in meta_cols]
    print('FRAME COLS:', frame_cols)
    
    y_pos = np.arange(len(plot_df))
    fig, ax = plt.subplots(figsize=(12, max(6, len(plot_df) * 0.4)))

    # Safe color mapping - use .map() which handles missing keys by returning NaN
    colors = plot_df["res_class"].map(COLOR_MAP)
    # Replace any NaN colors with a default color (gray)
    colors = colors.fillna('#95a5a6')
    
    ax.barh(y_pos, plot_df["mean"], color=colors,
            edgecolor="black", height=bar_height)

    if 'std' in plot_df.columns:
        ax.errorbar(plot_df["mean"], y_pos, xerr=plot_df["ci_half_width"],
                    fmt='none', ecolor='black', elinewidth=1.5, capsize=3, capthick=1.5)

    # --- Per-frame dots ---
    #if frame_cols:
    #    for i, row in plot_df.iterrows():
    #        vals = row[frame_cols].dropna().values
    #        ax.scatter(vals, np.full(len(vals), i),
    #                   color='black', s=8, zorder=5, alpha=0.6, linewidths=0)

    y_labels = [f"$\\bf{{{l}}}$" if i in top_indices else l
                for i, l in enumerate(plot_df["res_label"])]
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=14)
    ax.tick_params(axis='x', labelsize=20)
    ax.set_xlabel(r'$\lambda_{\max,\ \mathrm{native}} - \lambda_{\max,\ q=0}\ \mathrm{(nm)}$', fontsize=30)
    if x_min is not None:
        ax.set_xlim(left=x_min)
    if x_max is not None:
        ax.set_xlim(right=x_max)
    ax.set_title(title, fontsize=16, fontweight='bold')
    ax.axvline(0, color='black', linewidth=1, alpha=0.5)
    ax.invert_yaxis()
    
    # Only create legend for classes that actually appear in the data
    existing_classes = plot_df["res_class"].unique()
    handles = [plt.Line2D([0], [0], color=COLOR_MAP.get(k, '#95a5a6'), lw=8) for k in existing_classes]
    ax.legend(handles, existing_classes, fontsize=12)

    plt.tight_layout()
    plt.savefig(f"{fname}_grouped_bar.png", dpi=300, bbox_inches='tight')
    plt.close()
    return plot_df


def colored_mae_frames(df,mae_list,output_path):
    #make_colorbar
    mae_structures=[]    
    #for structure in mae_list:
    return

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Visualize residue contributions from QM/MM calculations')
    parser.add_argument('-pp', '--param_folders', nargs='+', required=True,help='Parameter folder(s) containing .txt files (e.g., ./parameters_000001 or ./parameters_*)')
    parser.add_argument('-o', '--output', type=str, default='contributions',help='Output file prefix (default: contributions)')
    parser.add_argument('-t','--title', type=str,default=None)
    parser.add_argument('--asl',type=str, default=None,help='optionally filter specific residues to plot')
    parser.add_argument('--xmin', type=float, default=None, help='X-axis minimum value')
    parser.add_argument('--xmax', type=float, default=None, help='X-axis maximum value')
    
    args = parser.parse_args()
    
    

    # Expand glob patterns if any
    param_folders = []
    for folder in args.param_folders:
        if '*' in folder or '?' in folder:
            param_folders.extend(sorted(glob.glob(folder)))
        else:
            param_folders.append(folder)
    print(os.path.dirname(param_folders[0]))
    mae_path = glob.glob(os.path.join(os.path.dirname(param_folders[0]), "*.mae"))[0]

    print(f"Found {len(param_folders)} parameter folders:")
    #for f in param_folders
        #print(f"  {f}")
    
    # Read MAE file
    if not os.path.exists(mae_path):
        print(f"Warning: MAE file not found: {mae_path}")
        print("Residue codes will be 'UNK'")
        mae_st = None
    else:
        mae_st = structure.StructureReader.read(mae_path)
    
    # Create matrix
    matrix = create_frame_contrib_matrix(param_folders, mae_st)
    if args.asl:
        from schrodinger.structutils import analyze
        asl_atoms = analyze.evaluate_asl(mae_st, args.asl)
        valid_residues = {(mae_st.atom[i].molecule_number, mae_st.atom[i].resnum) for i in asl_atoms}
        matrix = matrix[matrix.index.droplevel('rescode').isin(valid_residues)]
    
    if matrix.empty:
        print("No data found!")
        sys.exit(1)

    # Plot all residues with x limits
    plot_bar_grouped(matrix, f'{args.output}_all_residues', 
                    title=f'All Residue Contributions {args.title}', 
                    top_n=0, x_min=args.xmin, x_max=args.xmax)
    
    # Plot top 25 with x limits
    top25_indices = matrix['mean'].abs().nlargest(25).index
    matrix_top25 = matrix.loc[top25_indices]
    plot_bar_grouped(matrix_top25, f'{args.output}_top25_residues_{args.title}', 
                    title=f'Top 25 Residue Contributions {args.title}', 
                    top_n=0, x_min=args.xmin, x_max=args.xmax)
    
    # Save data
    matrix.reset_index().to_csv(f'{args.output}.csv', index=False)
    print(f"Processed {len(matrix)} residues")
    print(f"Saved: {args.output}.csv")