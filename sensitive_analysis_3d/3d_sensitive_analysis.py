import csv
import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

# Set Times New Roman font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'

def format_depth(depth):
    """Format depth value, remove negative sign and format to 2 decimal places"""
    return f"{abs(float(depth)):.2f}"

def format_afr(afr):
    """Format afr value"""
    return f"{float(afr):.2f}"

def format_threshold(thr):
    """Format threshold value"""
    return f"{float(thr):.2f}"

# Read CSV file
data = []
with open('3d_minpoint_detailed_results.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

print(f"Total rows read: {len(data)}")

# Group by (Lx, afr, threshold, timestep)
# We'll organize as: for each L, show all combinations of (afr, threshold) with all timesteps
grouped_by_L = defaultdict(lambda: defaultdict(list))

for row in data:
    Lx = int(float(row['Lx']))
    afr = float(row['afr'])
    threshold = float(row['threshold'])
    timestep = int(float(row['timestep']))
    depth = abs(float(row['min_z_scaled']))  # Remove negative sign
    
    key = (afr, threshold)
    grouped_by_L[Lx][key].append({
        'timestep': timestep,
        'depth': depth,
        'min_x': float(row['min_x']),
        'min_y': float(row['min_y']),
        'num_min_points': int(float(row['num_min_points']))
    })

# Sort timesteps within each group
for Lx in grouped_by_L:
    for key in grouped_by_L[Lx]:
        grouped_by_L[Lx][key].sort(key=lambda x: x['timestep'])

# Get unique values for sorting
L_values = sorted(grouped_by_L.keys())
afr_values = sorted(set(float(row['afr']) for row in data))
threshold_values = sorted(set(float(row['threshold']) for row in data))

print(f"L values: {L_values}")
print(f"afr values: {afr_values}")
print(f"threshold values: {threshold_values}")

# Create separate table for each L value
for L_value in L_values:
    print(f"\nProcessing L={L_value}...")
    
    # Prepare table data
    table_data = []
    depth_values = []  # For potential color coding
    
    # Sort by (afr, threshold)
    param_keys = sorted(grouped_by_L[L_value].keys())
    
    for param_idx, (afr, threshold) in enumerate(param_keys):
        timestep_data = grouped_by_L[L_value][(afr, threshold)]
        
        for ts_idx, ts_data in enumerate(timestep_data):
            depth = ts_data['depth']
            depth_values.append(depth)
            
            if ts_idx == 0:
                # First row for this parameter combination
                table_row = [
                    str(L_value * 10),  # L multiply by 10
                    format_afr(afr),
                    format_threshold(threshold),
                    str(ts_data['timestep']),
                    format_depth(depth)
                ]
            else:
                # Continuation rows - empty L, afr, threshold
                table_row = [
                    '',
                    '',
                    '',
                    str(ts_data['timestep']),
                    format_depth(depth)
                ]
            
            table_data.append(table_row)
    
    # Column headers
    headers = ['L', r'$a_{fr}^2$', r'$u_c$', r'$t$', 'Depth']
    
    # Create figure - more compact layout
    fig_height = len(table_data) * 0.25 + 0.65
    fig, ax = plt.subplots(figsize=(7.5, fig_height))
    ax.axis('tight')
    ax.axis('off')
    
    # Create table
    table = ax.table(cellText=table_data, colLabels=headers, cellLoc='center', loc='center')
    
    # Style the table - LaTeX three-line style, very compact
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1, 1.15)  # More compact row height
    
    # Remove all cell borders first and set compact padding
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.set_edgecolor('white')
        cell.set_facecolor('white')
        cell.PAD = 0.02  # Reduced padding for tighter spacing
    
    # Header styling - bold, no background color
    for i in range(len(headers)):
        cell = table[(0, i)]
        cell.set_text_props(weight='bold', color='black', fontsize=11)
        cell.set_facecolor('white')
    
    # Add only the three main horizontal lines (LaTeX booktabs style)
    # Top line (above header)
    for j in range(len(headers)):
        cell = table[(0, j)]
        cell.set_linewidth(0)
        cell.visible_edges = 'T'
        cell.set_edgecolor('black')
        cell.set_linewidth(1.5)
    
    # Middle line (below header)
    for j in range(len(headers)):
        cell = table[(1, j)]
        cell.visible_edges = 'T'
        cell.set_edgecolor('black')
        cell.set_linewidth(1.0)
    
    # Bottom line (below last row)
    for j in range(len(headers)):
        cell = table[(len(table_data), j)]
        cell.visible_edges = 'B'
        cell.set_edgecolor('black')
        cell.set_linewidth(1.5)
    
    # Add subtle lines between parameter groups (every 6 rows for different afr/threshold combos)
    timesteps_per_group = len(grouped_by_L[L_value][param_keys[0]])
    for i in range(1, len(table_data) + 1):
        if i % timesteps_per_group == 1 and i > 1:
            for j in range(len(headers)):
                cell = table[(i, j)]
                cell.visible_edges = 'T'
                cell.set_edgecolor('gray')
                cell.set_linewidth(0.5)
    
    # Adjust column widths - more compact spacing
    for i in range(len(table_data) + 1):
        table[(i, 0)].set_width(0.10)  # L column
        table[(i, 1)].set_width(0.12)  # afr column
        table[(i, 2)].set_width(0.12)  # threshold column
        table[(i, 3)].set_width(0.15)  # timestep column
        table[(i, 4)].set_width(0.18)  # depth column
    
    plt.tight_layout()
    filename = f'3d_depth_analysis_{L_value * 10}.png'
    plt.savefig(filename, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.02)
    print(f"✓ Table for L={L_value * 10} saved as '{filename}'")
    
    plt.close()

print("\n" + "="*80)
print("Summary:")
print(f"- Total data points: {len(data)}")
print(f"- L values: {[L*10 for L in L_values]}")
print(f"- afr values: {afr_values}")
print(f"- threshold values: {threshold_values}")
print(f"- Tables generated: {len(L_values)} (L={', '.join([str(L*10) for L in L_values])})")
print(f"- Table style: LaTeX three-line (booktabs), compact")
print(f"- Font: Times New Roman")
print(f"- Depth values: Absolute values (negative signs removed)")
print("="*80)

