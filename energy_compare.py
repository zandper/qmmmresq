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
import textscrape


def concat_param(folder_path):
    all_text = ""
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, 'r', encoding='utf-8') as file:
                all_text += file.read() + "\n"  # Add newline to separate files
    return(all_text)


def plot(df):
    # Hover text
    print(df.columns)
    df['hover'] = df.apply(
        lambda row: "<br>".join([f"{col}: {row[col]}" for col in df.columns]),
        axis=1
    )
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
    fig.write_html('bar.html')

    # Scatter plot using Plotly
    df.sort_values(by='cr_dist',ascending=True, inplace=True) # needs to be sorted for correct lines
    fig2 = px.scatter(
        df,
        x='cr_dist',            # x-axis: sidechain distance
        y='delta_lambda',       # y-axis: delta lambda
        labels={
            'cr_dist': 'Minimum Chromophore-Residue Distance (Å)',
            'delta_lambda': 'Δλ (nm)'
        },
        title='Delta Lambda vs. Minimum Chromophore-Residue Distance',
    )
    fig2.add_trace(
    go.Scatter(
        x=df['cr_dist'],
        y=df['delta_lambda'],
        mode='lines+markers',
        line_shape='linear',
        hovertext=df['hover'],
        hovertemplate="%{hovertext}<extra></extra>"
    ))

    # Save plot as HTML
    fig2.write_html('scatter.html')


#original energy
original_log_path = "./qsite_OCPO-pH-7.4.out"
o_energy = output.QSiteOutput(original_log_path).energy
o_lamb =textscrape.extract_first_wavelength(original_log_path)

folder_path = './parameters/'

# Convert text to DataFrame
df = pandas.read_csv(StringIO(concat_param(folder_path)), sep="\t")
df.columns = ['resnum', 'energy', 'lambda']
df = df.sort_values(by='resnum') #sort by residue number
df['delta_lambda'] = o_lamb - df['lambda'] # Compute Delta Lambda

if __name__ == "__main__":

    #Add other parameters
    mae_path = "/home/zap22001/charge_diff/qsite_OCPO-pH-7.4.mae"
    mae_st = structure.StructureReader.read(mae_path)

    mae_df=[]
    chrom_st = mae_st.molecule[2].extractStructure()
    #print(chrom_st)
    for res in mae_st.residue:
        #st1 is chromphore
        cr_dist = measure.get_shortest_distance(res.extractStructure(), st2=chrom_st)
        print(cr_dist[0])
        mae_df.append({    
            'resnum':res.resnum,
            'rescode':res.pdbres,
            'cr_dist':cr_dist[0]
        })

    mae_df = pandas.DataFrame(mae_df)

    df = df.merge(mae_df, on='resnum', how='left')
    print(df)
    plot(df)
    #df['abs_delta_lambda']= df['delta_lambda'].abs()
    #df.sort_values(by='abs_delta_lambda',ascending=False, inplace=True)
    #print(df[['resnum','rescode','delta_lambda','cr_dist']][:50])
    #scatter_plot(df)