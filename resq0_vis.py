#!/usr/bin/env python3
import os
from io import StringIO
import pandas
import numpy
import re
import math
import matplotlib as plt
from schrodinger import structure
from schrodinger.application.qsite import output
from schrodinger.structutils import analyze, measure
import plotly.express as px
import plotly.graph_objects as go
import plotly.io
plotly.io.renderers.default = 'browser'
import utils.textscrape
import argparse

from schrodinger import structure
from schrodinger.structutils import color
import matplotlib.cm as cm
import matplotlib.colors as mcolors

def concat_param_folder(folder_path): 
    all_text = "" 
    for filename in os.listdir(folder_path): 
        if filename.endswith(".txt"): 
            file_path = os.path.join(folder_path, filename) 
            with open(file_path, 'r', encoding='utf-8') as file: 
                all_text += file.read() + "\n" # Add newline to separate files return(all_text)
    print(all_text)
    return (all_text)

def concat_multiple_param_paths(param_paths):
    """
    Concatenate .txt files from multiple param paths and average values.
    
    Args:
        param_paths: List of folder paths containing parameter .txt files
    
    Returns:
        DataFrame with averaged values across all frames
    """
    all_dfs = []
    
    for path in param_paths:
        path = os.path.realpath(path)
        if not os.path.exists(path):
            print(f"Warning: Path does not exist: {path}")
            continue
            
        all_text = concat_param_folder(path)
        if not all_text:
            print(f"Warning: No .txt files found in {path}")
            continue
            
        df = pandas.read_csv(StringIO(all_text), sep="\t", header=None)
        df.columns = ['molnum', 'resnum', 'energy', 'lambda']
        df['frame'] = os.path.basename(path)  # Add frame identifier
        all_dfs.append(df)
    
    if not all_dfs:
        raise ValueError("No valid parameter folders found")
    
    # Combine all dataframes
    combined_df = pandas.concat(all_dfs, ignore_index=True)
    
    # Group by molnum and resnum to compute averages and std deviations
    grouped = combined_df.groupby(['molnum', 'resnum']).agg({
        'energy': ['mean', 'std'],
        'lambda': ['mean', 'std'],
        'frame': 'count'  # number of frames
    }).reset_index()
    
    # Flatten column names
    grouped.columns = ['molnum', 'resnum', 'energy_mean', 'energy_std', 
                       'lambda_mean', 'lambda_std', 'n_frames']
    
    # Convert to int
    grouped['resnum'] = grouped['resnum'].astype(int)
    grouped['molnum'] = grouped['molnum'].astype(int)
    
    return grouped, combined_df

def plot_bar(df,fname):
    # Plotly bar chart
    fig = px.bar(
        df,
        x='resnum',  # or another x-axis column
        y='delta_lambda',  # or another y-axis column
        labels={'resnum': 'Residue Number', 'delta_lambda': 'Δλ (nm)'},
        title="Interactive Bar Chart",
    )

    # Set hover text
    fig.update_traces(
    hovertemplate="%{hovertext}<extra></extra>",  # use your custom hover text
    hovertext=df['hover']
    )
    bar_fname = f'{fname}_bar.html'
    print(f'saving bargraph as{bar_fname}')
    fig.write_html(bar_fname)

def plot_scatter(df,fname):

    # Automatically detect the first cr_dist_* column
    cr_dist_cols = [col for col in df.columns if col.startswith('cr_dist')]
    if not cr_dist_cols:
        raise ValueError("No 'cr_dist_*' columns found in the dataframe.")
    cr_dist_col = cr_dist_cols[0]
    print(f"Using '{cr_dist_col}' for scatter plot")

    # Scatter plot using Plotly
    df.sort_values(by=cr_dist_col,ascending=True, inplace=True)
    fig2 = px.scatter(
        df,
        x=cr_dist_col,           
        y='delta_lambda', 
        labels={
            'cr_dist': 'Minimum Chromophore-Residue Distance (Å)',
            'delta_lambda': '$\Delta \lambda_{\mathrm{max}}$ (nm)'
        },
        title='Delta Lambda vs. Minimum Chromophore-Residue Distance',
    )
    fig2.add_trace(
    go.Scatter(
        x=df[cr_dist_col],
        y=df['delta_lambda'],
        mode='lines+markers',
        line_shape='linear',
        hovertext=df['hover'],
        hovertemplate="%{hovertext}<extra></extra>"
    ))
    # Save plot as HTML
    scatter_fname = f'{fname}_scatter.html'
    print(f'saving bargraph as{scatter_fname}')
    fig2.write_html(scatter_fname)

import matplotlib.pyplot as plt
import numpy as np


def plot_bar_grouped(df, fname, title, top_n=3, bar_height=0.8):
    """
    Horizontal bar plot colored by residue class.
    Highlights top N residues by absolute Δλ in bold on y-axis.
    """
    # Color map for residue classes
    color_map = {
        "(+) Charge": "#4a90e2",  # muted blue
        "(-) Charge": "#d64541",  # soft red
        "Polar": "#9b59b6",               # muted purple
        "Nonpolar": "#f1c40f",            # softer yellow
        "Other/Nonprotein": "#95a5a6"     # grayish soft
    }

    # Sort by class and chromophore distance
    df_sorted = df.sort_values(by=["res_class", "cr_dist"], ascending=[True, True]).reset_index(drop=True)

    # Identify top N residues by |delta_lambda|
    top_indices = df_sorted["delta_lambda"].abs().nlargest(top_n).index

    # Build plot with bigger figure and fonts
    fig, ax = plt.subplots(figsize=(12, max(6, len(df_sorted) * 0.4)))  # bigger width & height

    y_positions = np.arange(len(df_sorted))
    colors = df_sorted["res_class"].map(color_map)

    ax.barh(
        y_positions,
        df_sorted["delta_lambda"],
        color=colors,
        edgecolor="black",
        height=bar_height
    )

    # Y labels: bold top N
    y_labels = []
    for i, label in enumerate(df_sorted["res_label"]):
        if i in top_indices:
            y_labels.append(f"$\\bf{{{label}}}$")  # bold top N
        else:
            y_labels.append(label)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=1)  # bigger font

    # Axis labels & title
    ax.set_xlabel(r'$\lambda_{\mathrm{max, Native}} - \lambda_{\mathrm{max, q=0}}$ (nm)', fontsize=50)
    ax.set_ylabel("", fontsize=0)
    ax.set_title(title, fontsize=0, fontweight='bold')

    ax.invert_yaxis()  # closest residues on top
    ax.set_ylim(len(df)-0.4, -0.6)

    # Legend: bigger text
    handles = [
        plt.Line2D([0], [0], color=color_map[k], lw=8)
        for k in color_map if k in df_sorted["res_class"].unique()
    ]
    labels = [k for k in color_map if k in df_sorted["res_class"].unique()]
    legend = ax.legend(handles, labels, title="", fontsize=20,loc="upper right")
    legend.get_frame().set_linewidth(1.5)

    # Tick params
    ax.tick_params(axis='x', labelsize=30, width=1.5, length=6)
    ax.tick_params(axis='y', labelsize=25, width=1.5, length=6)

    plt.tight_layout()
    outname = f"{fname}_grouped_bar.png"
    print(f"saving matplotlib grouped bar as {outname}")
    plt.savefig(outname, dpi=300, bbox_inches='tight')
    plt.close()


def apply_lambda_coloring(mae_file, df):
    """
    Color atoms in the MAE file based on delta_lambda values in df.
    Iterates through the DataFrame by resnum and colors all atoms in that residue.
    Prints residue number, insertion code, and atoms colored.
    """
    st = next(structure.StructureReader(mae_file))

    # Blue-white-red colormap centered at 0, -15 blue, 0 white, 15 red
    max_abs = math.ceil(df["delta_lambda"].abs().max())
    if max_abs == 0:
        max_abs = 1.0

    norm = mcolors.Normalize(vmin=-max_abs, vmax=max_abs, clip=True)
    cmap = cm.get_cmap("bwr")

    colored_atoms = 0
    colored_residues = 0

    for _, row in df.iterrows():
        molnum = int(row["molnum"])
        df_resnum = int(row["resnum"])
        lam = float(row["delta_lambda"])

        r, g, b, _ = cmap(norm(lam))
        rgb255 = (int(r*255), int(g*255), int(b*255))

        # Find residue
        res_obj = next(
            (r for r in st.residue
             if r.molecule_number == molnum and int(str(r.resnum).strip()) == df_resnum),
            None
        )

        if res_obj is None:
            print(f"WARNING: Residue mol {molnum} res {df_resnum} not found")
            continue

        colored_residues += 1
        print(f"Coloring residue: mol {molnum} res {res_obj.resnum}")

        for atom in res_obj.atom:
            atom.color = rgb255
            colored_atoms += 1

    print(f"Total colored residues: {colored_residues}, total colored atoms: {colored_atoms}")
    colored_mae_filename="colored.mae" 
    st.write(colored_mae_filename) #("mae_file")
    print(f"Saved colored MAE: {colored_mae_filename}")

    # ---- Save colorbar with bigger/more legible font ----
    fig, ax = plt.subplots(figsize=(6, 0.5)) 
    fig.subplots_adjust(bottom=0.5)

    cb1 = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                       cax=ax, orientation='horizontal')

    # Increase label size
    cb1.set_label(r'$\lambda_{\mathrm{max, Native}} - \lambda_{\mathrm{max, q=0}}$ (nm)', fontsize=25)
    # Increase tick label size
    cb1.ax.tick_params(labelsize=14, width=2, length=8)
    # Optional: bold ticks
    for tick in cb1.ax.get_xticklabels():
        tick.set_fontweight('bold')

    plt.savefig("colorbar.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved colorbar.png")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-pp', '--param_paths', type=str, default='./parameters/', help='Path to folder of measured params')
    parser.add_argument('-a','--asl',type=str,default=None,help='Manually define residues in QM')
    parser.add_argument('-t','--title',type=str,default='')
    args = parser.parse_args()

    #  establish paths
    in_path = os.path.realpath(args.in_path)
    param_path = os.path.realpath(args.param_path)
    with open(in_path) as f:
        text = f.read()
    mae_path = (re.findall(r"MAEFILE:\s*(\S+)", text))[0]
    mae_fname = os.path.basename(mae_path)
    native_out_path = in_path.replace(".in",".out")

    # Get native energy / wavelength from native.log path
    o_energy = output.QSiteOutput(native_out_path).energy
    o_lamb = utils.textscrape.extract_first_wavelength(native_out_path)
    print(o_lamb)

    # Convert text to DataFrame
    alltext=folder(os.path.realpath(args.param_path))
    df = pandas.read_csv(StringIO(alltext), sep="\t", header=None)
    df.columns = ['molnum','resnum', 'energy', 'lambda']
    df['resnum'] = df['resnum'].astype(int)
    df['molnum'] = df['molnum'].astype(int)
    df = df.sort_values(by='resnum') #sort by residue number
    df['delta_lambda'] = o_lamb - df['lambda'] # Compute Delta Lambda
    mae_st = structure.StructureReader.read(mae_path)

    # ASL defining QM / chromophore
    if args.asl is None:
        raise ValueError("You must provide an ASL selection via --asl for the QM/chromophore region")

    chrom_indices = analyze.evaluate_asl(mae_st, args.asl)
    if not chrom_indices:
        raise RuntimeError(f"ASL '{args.asl}' returned zero atoms")
    chrom_st = mae_st.extract(chrom_indices)

    # Build distance dataframe using residues from parameter files
    mae_rows = []

    for _, row in df.iterrows():

        molnum = int(row["molnum"])
        resid  = int(row["resnum"])

        # Select residue by molecule number AND residue number
        res_asl = f"mol.num {molnum} AND res.num {resid}"
        res_indices = analyze.evaluate_asl(mae_st, res_asl)

        if not res_indices:
            raise RuntimeError(f"ASL '{res_asl}' returned zero atoms")

        res_st = mae_st.extract(res_indices)

        # shortest distance residue ↔ QM region
        cr_dist = measure.get_shortest_distance(res_st, st2=chrom_st)[0]

        # get residue identity safely
        res_obj = next(
            (r for r in mae_st.residue
            if r.molecule_number == molnum and r.resnum == resid),
            None
        )

        rescode = res_obj.pdbres if res_obj else "UNK"

        mae_rows.append({
            "molnum": molnum,
            "resnum": resid,
            "rescode": rescode.strip(),
            "cr_dist": cr_dist
        })

    mae_df = pandas.DataFrame(mae_rows)
    # merge with lambda dataframe
    df = df.merge(mae_df, on=["molnum", "resnum"], how="left")

    fname = in_path.replace(".in","")
    # Hover text
    df['hover'] = df.apply(lambda row: "<br>".join([f"{col}: {row[col]}" for col in df.columns]), axis=1)

    plot_bar(df, fname)
    plot_scatter(df, fname)
    # --- Residue classification ---
    pos_charged = {"ARG", "LYS", "HIS"}
    neg_charged = {"ASP", "GLU","ASH"}
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

    df["res_class"] = df["rescode"].apply(classify_residue)
    #df["res_label"] = ("Mol" + df["molnum"].astype(str) + "_" + df["rescode"] + df["resnum"].astype(str))
    df["res_label"] = df["rescode"] + df["resnum"].astype(str)
    class_order = [
    "(+) Charge",
    "(-) Charge",
    "Polar",
    "Nonpolar",
    "Other/Nonprotein"
    ]

    df["res_class"] = pandas.Categorical(
        df["res_class"],
        categories=class_order,
        ordered=True
    )
    df_sorted = df.sort_values(by='resnum')
    # Save to CSV without the index
    df_sorted.to_csv("debug.csv", index=False)
    apply_lambda_coloring(mae_path,df)
    df = df.sort_values(by=["res_class", "cr_dist"], ascending=[True, True])
    print(len(df_sorted))
    plot_bar_grouped(df,(fname+'_all'),args.title)
    df_top = df.reindex(df['delta_lambda'].abs().sort_values(ascending=False).index)[:25].copy()
    plot_bar_grouped(df_top,fname,args.title)
    #print(df.columns)
    #apply_lambda_coloring(args.mae,df)