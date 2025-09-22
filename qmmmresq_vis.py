#!/usr/bin/env python3
import os
from io import StringIO
import pandas
import numpy
import re
import matplotlib as plt
from schrodinger import structure
from schrodinger.application.qsite import output
from schrodinger.structutils import measure
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
                all_text += file.read() + "\n"  # Add newline to separate files
    return(all_text)


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
    #parser.add_argument('-n', '--nlog_path', type=str, required=True, help='Path to natively_charged_protein .log file')
    #parser.add_argument('-i', '--in_path', type=str, required=True, help='Path to .in file')
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
    df = pandas.read_csv(StringIO(concat_param(os.path.realpath(args.param_path))), sep="\t")
    df.columns = ['resnum', 'energy', 'lambda']
    df = df.sort_values(by='resnum') #sort by residue number
    df['delta_lambda'] = o_lamb - df['lambda'] # Compute Delta Lambda

    # Calculate distances from qm molecule
    mae_st = structure.StructureReader.read(mae_path)
    molids=textscrape.get_molids(in_path)
    
    mae_df={}
    # First pass: populate residue identity (assumes residue list is same across molids)
    residue_list = mae_st.residue
    mae_df['resnum'] = [res.resnum for res in residue_list]
    mae_df['rescode'] = [res.pdbres for res in residue_list]
    print(molids)
    # For each molid, calculate distances and add to data dictionary
    for molid in molids:
        chrom_st = mae_st.molecule[molid].extractStructure()
        dist_list = []
        
        for res in residue_list:
            cr_dist = measure.get_shortest_distance(res.extractStructure(), st2=chrom_st)
            dist_list.append(cr_dist[0])
        
        # Add this molid's distances as a new column
        mae_df[f'cr_dist_id_{molid}'] = dist_list

    # Create DataFrame from the dictionar

    mae_df = pandas.DataFrame(mae_df)
    print(mae_df)
    df = df.merge(mae_df, on='resnum', how='left')
    print(df)
    fname = in_path.replace(".in","")

    # Hover text
    df['hover'] = df.apply(
        lambda row: "<br>".join([f"{col}: {row[col]}" for col in df.columns]),
        axis=1
    )

    plot_bar(df,fname)
    plot_scatter(df,fname)