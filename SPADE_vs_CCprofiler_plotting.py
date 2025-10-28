import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


SPADE_R1_10min = pd.read_csv('SPADE_R1_10min_protein_features.csv')
SPADE_R2_10min = pd.read_csv('SPADE_proteinfeatures_R2_10min.csv')
SPADE_R3_10min = pd.read_csv('SPADE_proteinfeatures_R3_10min.csv')
SPADE_R1_30min = pd.read_csv('SPADE_proteinfeatures_R1_30.csv')
SPADE_R2_30min = pd.read_csv('SPADE_proteinfeatures_R2_30min.csv')
SPADE_R3_30min = pd.read_csv('SPADE_proteinfeatures_R3_30min.csv')

CCprofiler_R1_10min = pd.read_csv('CCprofiler_clean_R1_10min.csv')
CCprofiler_R2_10min = pd.read_csv('CCprofiler_clean_R2_10min.csv')
CCprofiler_R3_10min = pd.read_csv('CCprofiler_clean_R3_10min.csv')
CCprofiler_R1_30min = pd.read_csv('CCprofiler_clean_R1_30min.csv')
CCprofiler_R2_30min = pd.read_csv('CCprofiler_clean_R2_30min.csv')
CCprofiler_R3_30min = pd.read_csv('CCprofiler_clean_R3_30min.csv')

import re
patron= re.compile(r"^SPADE_R\d_\w+$")
SPADE_peak_counts = []
for var_name, var_value in globals().items():
    if patron.match(var_name):
        if hasattr(var_value, 'shape'):
           print(var_name)
           row_df = var_value.shape[0]
           SPADE_peak_counts.append(row_df)


import re
patron= re.compile(r"^CCprofiler_R\d_\w+$")
CCprofiler_peak_counts = []
for var_name, var_value in globals().items():
    if patron.match(var_name):
        if hasattr(var_value, 'shape'):
           print(var_name)
           row_df = var_value.shape[0]
           CCprofiler_peak_counts.append(row_df)

# number of peaks identified

labels = ['R1_10min', 'R2_10min', 'R3_10min', 'R1_30min', 'R2_30min', 'R3_30min']

# Create an array of x-axis positions for the bars
x = np.arange(len(labels))

# Set the width of the bars
width = 0.35
8
# Create the figure and axes
fig, ax = plt.subplots()

# Plot the first set of bars
rects1 = ax.bar(x - width/2, SPADE_peak_counts, width, label='SPADE')

# Plot the second set of bars
rects2 = ax.bar(x + width/2, CCprofiler_peak_counts, width, label='CCprofiler')

# Add some labels, a title, and a legend
ax.set_ylabel('Number of Rows')
ax.set_title('Comparison of Identified Peaks')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.legend()

# Display the plot
plt.show()

#Multiplicity

datasets = {
    ("SPADE", "10min", "R1"): SPADE_R1_10min,
    ("SPADE", "10min", "R2"): SPADE_R2_10min,
    ("SPADE", "10min", "R3"): SPADE_R3_10min,
    ("SPADE", "30min", "R1"): SPADE_R1_30min,
    ("SPADE", "30min", "R2"): SPADE_R2_30min,
    ("SPADE", "30min", "R3"): SPADE_R3_30min,
    ("CCprofiler", "10min", "R1"): CCprofiler_R1_10min,
    ("CCprofiler", "10min", "R2"): CCprofiler_R2_10min,
    ("CCprofiler", "10min", "R3"): CCprofiler_R3_10min,
    ("CCprofiler", "30min", "R1"): CCprofiler_R1_30min,
    ("CCprofiler", "30min", "R2"): CCprofiler_R2_30min,
    ("CCprofiler", "30min", "R3"): CCprofiler_R3_30min,
}

df = pd.concat(datasets, axis=1)


protdf = df.xs('protein_id', level= 3, axis=1)
df_long = protdf.stack(level=[0, 1,2])

# Rename the index levels for clarity
df_long.index.names = ['original_index', 'group', 'time', 'replicate']

# The values of the DataFrame (protein IDs) become a Series,
# so we convert it back to a DataFrame and name the column 'protein_id'
df_long = df_long.reset_index(name='protein_id')

# Count the occurrences of each unique protein ID for each time and replicate
counts = df_long.groupby(['group','time', 'replicate', 'protein_id']).size().reset_index(name='count')

print("Protein ID counts per time and replicate:")
print(counts)

plt.figure(figsize=(10, 6))

# Plot the histograms
sns.histplot(data=counts, x='count', hue='group', multiple='dodge', bins=20, shrink=0.8)

plt.title('Protein Multiplicity')
plt.xlabel('Number of Occurrences')
plt.ylabel('Number of Protein IDs')
plt.show()

#BOXPLOT OF APEX 

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannotations.Annotator import Annotator
df_long = apex.stack(level=[0, 1,2,3])
# Rename the index levels for clarity
df_long.index.names = ['original_index', 'group', 'time', 'replicate','protein_id','fraction','intensity']
df_long

df_long = apex.stack(level=[0, 1,2,3])
# Rename the index levels for clarity
df_long.index.names = ['original_index', 'group', 'time', 'replicate','apex_location']
df_long = df_long.reset_index(name='apex_value')
df_long_clean = df_long.loc[:,('group','time','replicate','apex_value')]
df_long_clean

time1_data = df_long_clean[df_long_clean['group'] == 'CCprofiler']['apex_value']
time2_data = df_long_clean[df_long_clean['group'] == 'SPADE']['apex_value']
statistic, p_value = stats.mannwhitneyu(time1_data, time2_data)
#statistic, p_value = stats.ttest_ind(time1_data, time2_data, equal_var=False)
print("\nIndependent T-Test Results (Welch's):")
print(f"Statistic: {statistic}")
print(f"P-value: {p_value}")
sns.set_style("whitegrid")

# Define the pairs to compare
pairs = [ (('10min', 'CCprofiler'), ('10min', 'SPADE')),(('30min', 'CCprofiler'), ('30min', 'SPADE'))]

# Create the boxplot
plt.figure(figsize=(10, 6))
ax = sns.boxplot(
    data=df_long_clean,
    x='time',
    y='apex_value',
    hue='group',
    width=0.6  )


# Use the Annotator to add the statistical significance
annotator = Annotator(
    ax,
    pairs,
    data=df_long_clean,
    x='time',
    y='apex_value',
    hue= 'group'
)

# Add the p-value with the Mann-Whitney U test
# You can customize the p-value format using the pvalue_format argument
annotator.configure(test='Mann-Whitney', text_format='star', loc='outside')
annotator.apply_and_annotate()

# Add plot titles and labels
plt.xlabel('Apex Value Distribution', fontsize=12)
plt.ylabel('Apex Value', fontsize=12)

# Show the legend
plt.legend(title='Experimental Condition', loc='upper right')

# Display the plot
plt.show()

#lines of peaks 

df.loc[:, (slice(None),slice(None) , ['R3'], ['protein_id','fraction','intensity'])]

# Define list of protein IDs you want to keep
proteins_of_interest = ['P60710', 'P17182', 'P10649', 'Q9DBJ1', 'P08113']

# Extract protein_id column across all dataset/time combinations
protein_ids = peaks.xs('protein_id', level=3, axis=1)
mask = protein_ids.apply(lambda col: col.isin(proteins_of_interest))
filtered_peaks = peaks[mask.any(axis=1)]
filtered_peaks


##MATCHING IN MULTIPLICITY
pivoted_df = counts.pivot_table(index='protein_id', columns='group', values='count', fill_value=0).reset_index()

# Rename columns for clarity
pivoted_df.columns.name = None
pivoted_df = pivoted_df.rename(columns={'CCprofiler': 'CCprofiler_count', 'SPADE': 'SPADE_count'})


# 2. Identify proteins with matching counts
possible_same = (
    pivoted_df[pivoted_df['CCprofiler_count'] == pivoted_df['SPADE_count']]
    ['protein_id']
    .tolist()
)

# 3. Identify proteins with differing counts and determine the highest group
unique_or_different_df = pivoted_df[pivoted_df['CCprofiler_count'] != pivoted_df['SPADE_count']].copy()

# Create a new column to hold the name of the group with the highest count
#unique_or_different_df['highest_count_group'] = unique_or_different_df.apply(
#    lambda row: 'CCprofiler' if row['CCprofiler_count'] > row['SPADE_count'] else 'SPADE',
#    axis=1
#)

# 4. Save the results to the 'unique' vector
#unique = list(unique_or_different_df['protein_id']
                  #, unique_or_different_df['highest_count_group']))

print("Possible Same (Matching IDs and Counts):")
print(possible_same)

#plotting peaks

CCprofiler_R1_10min = pd.read_csv('R1_10min_prot.csv')
CCprofiler_R1_30min = pd.read_csv('R1_30min_prot.csv')
CCprofiler_R2_10min = pd.read_csv('R2_10min_prot.csv')
CCprofiler_R2_30min = pd.read_csv('R2_30min_prot.csv')
CCprofiler_R3_10min = pd.read_csv('R3_10min_prot.csv')
CCprofiler_R3_30min = pd.read_csv('R3_30min_prot.csv')
SPADE_R1_10min= pd.read_csv('SPADE_smointensity_R1_10.csv')
SPADE_R2_10min= pd.read_csv('SPADE_smointensity_R2_10.csv')
SPADE_R3_10min= pd.read_csv('SPADE_smointensity_R3_10.csv')
SPADE_R1_30min= pd.read_csv('SPADE_smointensity_R1_30.csv')
SPADE_R2_30min= pd.read_csv('SPADE_smointensity_R2_30.csv')
SPADE_R3_30min= pd.read_csv('SPADE_smointensity_R3_10.csv')


a = CCprofiler_R1_10min.nlargest(n= 81831, columns= 'intensity')
j= a['protein_id'].drop_duplicates()
top10= [j.head(n= 10)]
top10

##HEATMAP IDEA


# 2. Extract the data, labels, and ticks for the plot
matrix = heatmap_data
protein_ids = heatmap_data.index
groups = heatmap_data.columns.to_list()

# 3. Create the imshow plot
fig, ax = plt.subplots(figsize=(8, 10))
im = ax.imshow(matrix, cmap='YlGnBu')  # Use a colormap like 'YlGnBu' or 'viridis'

# 4. Set up the ticks and labels for the axes
ax.set_xticks(np.arange(len(groups)))
ax.set_yticks(np.arange(len(protein_ids)))
ax.set_xticklabels(groups)
ax.set_yticklabels(protein_ids)

# 5. Add a colorbar to act as a legend for the count values
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel("Count", rotation=-90, va="bottom")

# 6. Add a title
ax.set_title("Protein Counts by Group")

# 7. Ensure labels are correctly displayed
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

# 8. Display the plot
plt.tight_layout()
plt.show()

datasets = {
    ("SPADE", "10min", "R1"): SPADE_R1_10min,
    ("SPADE", "10min", "R2"): SPADE_R2_10min,
    ("SPADE", "10min", "R3"): SPADE_R3_10min,
    ("SPADE", "30min", "R1"): SPADE_R1_30min,
    ("SPADE", "30min", "R2"): SPADE_R2_30min,
    ("SPADE", "30min", "R3"): SPADE_R3_30min,
    ("CCprofiler", "10min", "R1"): CCprofiler_R1_10min,
    ("CCprofiler", "10min", "R2"): CCprofiler_R2_10min,
    ("CCprofiler", "10min", "R3"): CCprofiler_R3_10min,
    ("CCprofiler", "30min", "R1"): CCprofiler_R1_30min,
    ("CCprofiler", "30min", "R2"): CCprofiler_R2_30min,
    ("CCprofiler", "30min", "R3"): CCprofiler_R3_30min,
}

df = pd.concat(datasets, axis=1)

peaks = df.loc[:, (slice(None),slice(None) , ['R3'], ['protein_id','fraction','intensity'])]

# Define list of protein IDs you want to keep
proteins_of_interest = ['P60710', 'P17182', 'P10649', 'Q9DBJ1', 'P08113']

# Extract protein_id column across all dataset/time combinations
protein_ids = peaks.xs('protein_id', level=3, axis=1)
mask = protein_ids.apply(lambda col: col.isin(proteins_of_interest))
filtered_peaks = peaks[mask.any(axis=1)]
filtered_peaks

# 3. Create the line plot
plt.figure(figsize=(12, 8))
plt.title("Protein Intensity vs. Fraction")
plt.xlabel("Fraction")
plt.ylabel("Intensity")

# Get unique combinations of protein_id and method to plot each line
unique_lines = groupby(['protein_id', 'method'])
for (protein_id, method), group in unique_lines:
    label = f"{protein_id} ({method})"
    plt.plot(group['fraction'], group['intensity'], marker='o', label=label)

plt.legend(title='Protein and Method', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()