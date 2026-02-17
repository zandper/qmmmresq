#!/usr/bin/env python3
import os
from io import StringIO
import pandas
import numpy
import re
import matplotlib as plt
from schrodinger import structure
from schrodinger.application.qsite import output
from schrodinger.structutils import analyze, measure
import plotly.express as px
import plotly.graph_objects as go
import plotly.io
plotly.io.renderers.default = 'browser'
from utils import textscrape
import argparse

def concat_param(folder_path): 
    all_text = "" 
    for filename in os.listdir(folder_path): 
        if filename.endswith(".txt"): 
            file_path = os.path.join(folder_path, filename) 
            with open(file_path, 'r', encoding='utf-8') as file: 
                all_text += file.read() + "\n" # Add newline to separate files return(all_text)
    return (all_text)

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
            'delta_lambda': 'Δλ (nm)'
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

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Process MAE and input files")
    parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
    parser.add_argument('-pp', '--param_path', type=str, default='./parameters/', help='Path to folder of measured params')
    parser.add_argument('-a','--asl',type=str,default=None,help='Manually define residues in QM')
    args = parser.parse_args()

    #  establish paths
    in_path = os.path.realpath(args.in_path)
    param_path = os.path.realpath(args.param_path)
    mae_path = in_path.replace(".in",".mae")
    native_out_path = in_path.replace(".in",".out")

    # Get native energy / wavelength from native.log path
    o_energy = output.QSiteOutput(native_out_path).energy
    o_lamb =textscrape.extract_first_wavelength(native_out_path)
    print(o_lamb)

    # Convert text to DataFrame
    alltext= concat_param(os.path.realpath(args.param_path))
    df = pandas.read_csv(StringIO(alltext), sep="\t")
    df.columns = ['resnum', 'energy', 'lambda']
    df['resnum'] = df['resnum'].astype(int)
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

    for resid in df['resnum'].unique():

        # residue ASL
        print('RESIDUE ID',resid)
        res_asl = f"res.num {int(resid)}"

        res_indices = analyze.evaluate_asl(mae_st, res_asl)
        if not res_indices:
            raise RuntimeError(f"ASL '{res_asl}' returned zero atoms")
        res_st = mae_st.extract(res_indices)
                                      

        # shortest distance residue ↔ QM region
        cr_dist = measure.get_shortest_distance(res_st, st2=chrom_st)[0]

        # get residue identity
        res_obj = next((r for r in mae_st.residue if r.resnum == resid), None)
        rescode = res_obj.pdbres if res_obj else "UNK"

        mae_rows.append({
            'resnum': int(resid),
            'rescode': rescode,
            'cr_dist': cr_dist
        })

    mae_df = pandas.DataFrame(mae_rows)

    print(mae_df)

    # merge with lambda dataframe
    df = df.merge(mae_df, on='resnum', how='left')

    print(df)
    fname = in_path.replace(".in","")

    # Hover text
    df['hover'] = df.apply(
        lambda row: "<br>".join([f"{col}: {row[col]}" for col in df.columns]),
        axis=1
    )

    plot_bar(df, fname)
    plot_scatter(df, fname)
