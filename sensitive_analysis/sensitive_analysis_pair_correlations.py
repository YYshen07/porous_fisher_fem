import csv
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

# Set Times New Roman font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'

def format_pvalue(p):
    """Format p-value in scientific notation if very small"""
    p = float(p)
    if p < 0.001:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"

def format_correlation(r):
    """Format correlation coefficient to 3 decimal places"""
    return f"{float(r):.3f}"

def pair_type_name(pt, start, end):
    """Convert pair_type to readable format with step range"""
    start = int(float(start))
    end = int(float(end))
    mapping = {
        'minus20': f'$T_{{50}}-20$ ({start}-{end})',
        'plus20': f'$T_{{50}}+20$ ({start}-{end})',
        'minus_2pct_full': f'$T_{{50}}-0.02T_{{full}}$ ({start}-{end})',
        'plus_2pct_full': f'$T_{{50}}+0.02T_{{full}}$ ({start}-{end})'
    }
    return mapping.get(pt, pt)

def calculate_T50(rows_sorted):
    """Calculate T50 from the sorted rows"""
    # T50 is the end of minus20 or start of plus20 (they should be the same)
    for row in rows_sorted:
        if row['pair_type'] == 'minus20':
            return int(float(row['step_end']))
        elif row['pair_type'] == 'plus20':
            return int(float(row['step_start']))
    return None

# Read CSV file
data = []
with open('sensitive_analysis_pair_correlations.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        data.append(row)

# Group by (Lx, afr, threshold)
grouped = defaultdict(list)
for row in data:
    key = (int(float(row['Lx'])), float(row['afr']), float(row['threshold']))
    grouped[key].append(row)

# Sort groups
sorted_groups = sorted(grouped.items(), key=lambda x: (x[0][0], x[0][1], x[0][2]))

# Prepare table data - compact format with p-values on same line
table_data = []
pearson_values = []  # Store pearson r values for color coding
spearman_values = []  # Store spearman rho values for color coding
order = ['minus20', 'plus20', 'minus_2pct_full', 'plus_2pct_full']

for params, rows in sorted_groups:
    Lx, afr, threshold = params
    # Sort rows by pair_type
    rows_sorted = sorted(rows, key=lambda x: order.index(x['pair_type']))
    
    # Calculate T50 for this parameter set
    T50 = calculate_T50(rows_sorted)
    
    for idx, row in enumerate(rows_sorted):
        pearson_r = float(row['pearson_r'])
        spearman_r = float(row['spearman_r'])
        
        # Combine correlation and p-value on same line
        pearson_str = f"{format_correlation(row['pearson_r'])} ({format_pvalue(row['pearson_p'])})"
        spearman_str = f"{format_correlation(row['spearman_r'])} ({format_pvalue(row['spearman_p'])})"
        
        if idx == 0:
            table_row = [
                str(Lx * 10),  # L multiply by 10
                f"{afr:.2f}",
                f"{threshold:.2f}",
                str(T50),  # T50 column
                pair_type_name(row['pair_type'], row['step_start'], row['step_end']),
                pearson_str,
                spearman_str
            ]
        else:
            table_row = [
                '',
                '',
                '',
                '',  # Empty T50 for continuation rows
                pair_type_name(row['pair_type'], row['step_start'], row['step_end']),
                pearson_str,
                spearman_str
            ]
        table_data.append(table_row)
        pearson_values.append(pearson_r)
        spearman_values.append(spearman_r)

# Column headers - include L, T50 and afr with square
headers = ['L', r'$a^2$', r'$u_c$', r'$T_{50}$', 'Time Window', 'Pearson r (P value)', r'Spearman $\rho$ (P value)']

# Group data by L value
data_by_L = defaultdict(lambda: {'data': [], 'pearson': [], 'spearman': []})
idx = 0
for params, rows in sorted_groups:
    Lx, afr, threshold = params
    rows_sorted = sorted(rows, key=lambda x: order.index(x['pair_type']))
    
    for row_idx, row in enumerate(rows_sorted):
        L_key = Lx * 10  # 400, 600, 800
        if row_idx == 0:
            table_row = [
                str(L_key),  # L value
                f"{afr:.2f}",
                f"{threshold:.2f}",
                table_data[idx][3],  # T50 value
                pair_type_name(row['pair_type'], row['step_start'], row['step_end']),
                table_data[idx][5],  # Pearson string
                table_data[idx][6]   # Spearman string
            ]
        else:
            table_row = [
                '',  # Empty L for continuation rows
                '',
                '',
                '',  # Empty T50 for continuation rows
                pair_type_name(row['pair_type'], row['step_start'], row['step_end']),
                table_data[idx][5],  # Pearson string
                table_data[idx][6]   # Spearman string
            ]
        data_by_L[L_key]['data'].append(table_row)
        data_by_L[L_key]['pearson'].append(pearson_values[idx])
        data_by_L[L_key]['spearman'].append(spearman_values[idx])
        idx += 1

# Create separate table for each L value
total_red_pearson = 0
total_red_spearman = 0

for L_value in sorted(data_by_L.keys()):
    table_data_L = data_by_L[L_value]['data']
    pearson_values_L = data_by_L[L_value]['pearson']
    spearman_values_L = data_by_L[L_value]['spearman']
    
    # Create figure - very compact layout
    fig, ax = plt.subplots(figsize=(10, len(table_data_L) * 0.28 + 0.9))
    ax.axis('tight')
    ax.axis('off')
    
    # Add title with negative pad to bring it closer
    plt.title(f'Sensitivity Analysis for L={L_value} μm', fontsize=14, weight='bold', pad=-20, y=0.98)
    
    # Create table
    table = ax.table(cellText=table_data_L, colLabels=headers, cellLoc='center', loc='center')
    
    # Style the table - LaTeX three-line style, very compact
    table.auto_set_font_size(False)
    table.set_fontsize(10)  # Balanced font size
    table.scale(1, 1.2)  # Compact row height
    
    # Remove all cell borders first and set compact padding
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.set_edgecolor('white')
        cell.set_facecolor('white')
        cell.PAD = 0.03  # Reduce cell padding for more compact layout
    
    # Header styling - bold, no background color
    for i in range(len(headers)):
        cell = table[(0, i)]
        cell.set_text_props(weight='bold', color='black', fontsize=11)
        cell.set_facecolor('white')
    
    # Data cells - apply red color to values < 0.9
    for i in range(1, len(table_data_L) + 1):
        row_idx = i - 1
        for j in range(len(headers)):
            cell = table[(i, j)]
            
            # Check if this is a correlation column and value < 0.9
            if j == 5:  # Pearson column (now column index 5)
                if pearson_values_L[row_idx] < 0.9:
                    cell.set_text_props(color='red')
            elif j == 6:  # Spearman column (now column index 6)
                if spearman_values_L[row_idx] < 0.9:
                    cell.set_text_props(color='red')
    
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
        cell = table[(len(table_data_L), j)]
        cell.visible_edges = 'B'
        cell.set_edgecolor('black')
        cell.set_linewidth(1.5)
    
    # Add subtle lines between parameter groups (every 4 rows for different afr/u_c combos)
    for i in range(1, len(table_data_L) + 1):
        if i % 4 == 1 and i > 1:
            for j in range(len(headers)):
                cell = table[(i, j)]
                cell.visible_edges = 'T'
                cell.set_edgecolor('gray')
                cell.set_linewidth(0.5)
    
    # Adjust column widths - very compact spacing
    for i in range(len(table_data_L) + 1):
        table[(i, 0)].set_width(0.05)  # L column
        table[(i, 1)].set_width(0.06)  # afr column
        table[(i, 2)].set_width(0.06)  # u_c column
        table[(i, 3)].set_width(0.07)  # T50 column
        table[(i, 4)].set_width(0.28)  # Time Window column
        table[(i, 5)].set_width(0.24)  # Pearson column
        table[(i, 6)].set_width(0.24)  # Spearman column
    
    plt.tight_layout()
    filename = f'sensitive_analysis_pair_correlations_{L_value}.png'
    plt.savefig(filename, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.02)
    print(f"✓ Table for L={L_value} saved as '{filename}'")
    
    # Count red values
    red_p = sum(1 for v in pearson_values_L if v < 0.9)
    red_s = sum(1 for v in spearman_values_L if v < 0.9)
    total_red_pearson += red_p
    total_red_spearman += red_s
    
    plt.close()

print("\n" + "="*80)
print("Summary:")
print(f"- Total parameter combinations: {len(sorted_groups)}")
print(f"- Tables generated: 3 (L=400, 600, 800)")
print(f"- Pearson r < 0.9 (marked in red): {total_red_pearson}/{len(pearson_values)}")
print(f"- Spearman ρ < 0.9 (marked in red): {total_red_spearman}/{len(spearman_values)}")
print(f"- Table style: LaTeX three-line (booktabs), compact")
print(f"- Font: Times New Roman")
print("="*80)
