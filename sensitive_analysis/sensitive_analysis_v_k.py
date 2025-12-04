import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
from collections import defaultdict

# Set Times New Roman font
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'

def format_afr(value):
    """Format afr value to match filename convention"""
    txt = f"{value:.2f}"
    txt = txt.rstrip('0').rstrip('.')
    return txt

def format_threshold(value):
    """Format threshold value to match filename convention"""
    txt = f"{value:.2f}"
    txt = txt.rstrip('0').rstrip('.')
    return txt

def load_crack_data(Lx, Ly, afr, threshold, timestep):
    """Load crack tip data for a specific timestep"""
    afr_str = format_afr(afr)
    threshold_str = format_threshold(threshold)
    filename = f'Lx{Lx}_Ly{Ly}_afr{afr_str}_thr{threshold_str}_flag1_timestep{timestep}.csv'
    
    if not os.path.exists(filename):
        return None
    
    data = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            data.append({
                'x': float(row['x']),
                'y': float(row['y']),
                'angle': float(row['angle']),
                'curvature': float(row['curvature'])
            })
    
    # Sort by angle to ensure points are in consistent order
    data.sort(key=lambda p: p['angle'])
    
    return data

def calculate_point_velocities(data_prev, data_curr):
    """Calculate velocity for each point between two consecutive timesteps
    Returns arrays of velocities and corresponding curvatures
    Similar to collect_pair_data in the notebook"""
    if data_prev is None or data_curr is None:
        return None, None
    
    if len(data_prev) != len(data_curr):
        return None, None
    
    if len(data_prev) == 0:
        return None, None
    
    # Check if angles match (points should be in same order)
    angles_match = np.allclose([p['angle'] for p in data_prev], 
                               [p['angle'] for p in data_curr])
    if not angles_match:
        return None, None
    
    velocities = []
    curvatures = []
    
    # Calculate displacement for each corresponding point
    for p_prev, p_curr in zip(data_prev, data_curr):
        dx = p_curr['x'] - p_prev['x']
        dy = p_curr['y'] - p_prev['y']
        distance = np.sqrt(dx**2 + dy**2)
        
        # Use curvature from previous timestep (like in the notebook)
        k = abs(p_prev['curvature'])
        
        # Only include points with positive velocity and curvature
        if distance > 0 and k > 0:
            velocities.append(distance)
            curvatures.append(k)
    
    if len(velocities) > 0:
        return np.array(velocities), np.array(curvatures)
    return None, None

def calculate_vk_analysis(Lx, Ly, afr, threshold, step_start, step_end):
    """Calculate v-k analysis including correlations and slope for a time window
    SAME AS pair_correlations: calculate between TWO specific timesteps (start and end)
    Returns: correlations between curvature and distance, plus v-k slope"""
    
    # Load data for the start and end timesteps ONLY (like pair_correlations.py)
    data_start = load_crack_data(Lx, Ly, afr, threshold, step_start)
    data_end = load_crack_data(Lx, Ly, afr, threshold, step_end)
    
    if data_start is None:
        print(f"    Warning: No data for step {step_start}")
        return None
    
    if data_end is None:
        print(f"    Warning: No data for step {step_end}")
        return None
    
    # Calculate velocities for all points between start and end
    velocities, curvatures = calculate_point_velocities(data_start, data_end)
    
    if velocities is None or curvatures is None:
        print(f"    Warning: Failed to calculate velocities")
        return None
    
    if len(velocities) < 10:
        print(f"    Warning: Insufficient data points ({len(velocities)} collected)")
        return None
    
    # Calculate Pearson and Spearman correlations (curvature vs distance/velocity)
    # This should match pair_correlations.py exactly
    pearson_r, pearson_p = stats.pearsonr(curvatures, velocities)
    spearman_r, spearman_p = stats.spearmanr(curvatures, velocities)
    
    # 1. Linear fit: v = a·κ + b (direct linear relationship)
    linear_slope = None
    linear_r2 = None
    
    if len(velocities) > 10:
        slope, intercept, r_value, p_value, std_err = stats.linregress(curvatures, velocities)
        linear_slope = slope
        linear_r2 = r_value**2
    
    # 2. Power-law fit: v ∝ κ^n (via log-log regression)
    # Physical basis: Crack propagation speed often follows power-law with curvature
    # Taking log: log(v) = log(A) + n·log(κ)
    valid_mask = (velocities > 0) & (curvatures > 0)
    v_valid = velocities[valid_mask]
    k_valid = curvatures[valid_mask]
    
    power_exponent = None
    power_r2 = None
    
    if len(v_valid) > 10:
        log_v = np.log(v_valid)
        log_k = np.log(k_valid)
        
        # Linear regression in log-log space: log(v) = n·log(k) + log(A)
        # The slope 'n' is the power-law exponent
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_k, log_v)
        power_exponent = slope
        power_r2 = r_value**2
    
    return {
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'linear_slope': linear_slope,
        'linear_r2': linear_r2,
        'power_exponent': power_exponent,
        'power_r2': power_r2,
        'n_points': len(velocities)
    }

# Read summary file
summary_data = []
with open('sensitive_analysis_summary.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        summary_data.append(row)

# Calculate v-k slopes for each configuration
results = []

print("Processing configurations...")
for row in summary_data:
    Lx = int(float(row['Lx']))
    Ly = int(float(row['Ly']))
    afr = float(row['afr'])
    threshold = float(row['threshold'])
    fifty_step = int(float(row['fifty_step']))
    step_minus_2pct = int(float(row['step_minus_2pct_full']))
    step_plus_2pct = int(float(row['step_plus_2pct_full']))
    
    print(f"\nProcessing: L={Lx}, afr={afr}, u_c={threshold}, T50={fifty_step}")
    
    # Calculate analysis for T50 - 0.02%full to T50
    print(f"  Before T50 window: {step_minus_2pct} to {fifty_step} ({fifty_step - step_minus_2pct + 1} steps)")
    result_minus = calculate_vk_analysis(Lx, Ly, afr, threshold, step_minus_2pct, fifty_step)
    if result_minus:
        print(f"    → Pearson r: {result_minus['pearson_r']:.4f}, p: {result_minus['pearson_p']:.2e}")
        print(f"    → Spearman ρ: {result_minus['spearman_r']:.4f}, p: {result_minus['spearman_p']:.2e}")
        if result_minus['linear_slope'] is not None:
            print(f"    → Linear slope (v=a·κ+b): {result_minus['linear_slope']:.4f}, R²: {result_minus['linear_r2']:.4f}")
        if result_minus['power_exponent'] is not None:
            print(f"    → Power exponent n (v∝κ^n): {result_minus['power_exponent']:.4f}, R²: {result_minus['power_r2']:.4f}")
        print(f"    → Data points: {result_minus['n_points']}")
    else:
        print(f"    → Failed to calculate (insufficient data)")
    
    # Calculate analysis for T50 to T50 + 0.02%full
    print(f"  After T50 window: {fifty_step} to {step_plus_2pct} ({step_plus_2pct - fifty_step + 1} steps)")
    result_plus = calculate_vk_analysis(Lx, Ly, afr, threshold, fifty_step, step_plus_2pct)
    if result_plus:
        print(f"    → Pearson r: {result_plus['pearson_r']:.4f}, p: {result_plus['pearson_p']:.2e}")
        print(f"    → Spearman ρ: {result_plus['spearman_r']:.4f}, p: {result_plus['spearman_p']:.2e}")
        if result_plus['linear_slope'] is not None:
            print(f"    → Linear slope (v=a·κ+b): {result_plus['linear_slope']:.4f}, R²: {result_plus['linear_r2']:.4f}")
        if result_plus['power_exponent'] is not None:
            print(f"    → Power exponent n (v∝κ^n): {result_plus['power_exponent']:.4f}, R²: {result_plus['power_r2']:.4f}")
        print(f"    → Data points: {result_plus['n_points']}")
    else:
        print(f"    → Failed to calculate (insufficient data)")
    
    # Store results for both windows
    results.append({
        'Lx': Lx,
        'afr': afr,
        'threshold': threshold,
        'T50': fifty_step,
        'window_type': 'minus_2pct_full',
        'step_start': step_minus_2pct,
        'step_end': fifty_step,
        'pearson_r': result_minus['pearson_r'] if result_minus else None,
        'pearson_p': result_minus['pearson_p'] if result_minus else None,
        'spearman_r': result_minus['spearman_r'] if result_minus else None,
        'spearman_p': result_minus['spearman_p'] if result_minus else None,
        'linear_slope': result_minus['linear_slope'] if result_minus else None,
        'linear_r2': result_minus['linear_r2'] if result_minus else None,
        'power_exponent': result_minus['power_exponent'] if result_minus else None,
        'power_r2': result_minus['power_r2'] if result_minus else None,
        'n_points': result_minus['n_points'] if result_minus else None,
    })
    
    results.append({
        'Lx': Lx,
        'afr': afr,
        'threshold': threshold,
        'T50': fifty_step,
        'window_type': 'plus_2pct_full',
        'step_start': fifty_step,
        'step_end': step_plus_2pct,
        'pearson_r': result_plus['pearson_r'] if result_plus else None,
        'pearson_p': result_plus['pearson_p'] if result_plus else None,
        'spearman_r': result_plus['spearman_r'] if result_plus else None,
        'spearman_p': result_plus['spearman_p'] if result_plus else None,
        'linear_slope': result_plus['linear_slope'] if result_plus else None,
        'linear_r2': result_plus['linear_r2'] if result_plus else None,
        'power_exponent': result_plus['power_exponent'] if result_plus else None,
        'power_r2': result_plus['power_r2'] if result_plus else None,
        'n_points': result_plus['n_points'] if result_plus else None,
    })

# Save results to CSV
with open('sensitive_analysis_vk_slopes.csv', 'w', newline='') as f:
    fieldnames = ['Lx', 'afr', 'threshold', 'T50', 'window_type', 
                  'step_start', 'step_end',
                  'pearson_r', 'pearson_p', 'spearman_r', 'spearman_p',
                  'linear_slope', 'linear_r2', 
                  'power_exponent', 'power_r2', 
                  'n_points']
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(results)

print("\n" + "="*80)
print("Results saved to 'sensitive_analysis_vk_slopes.csv'")
print("="*80)

# Now create the table visualization
print("\nGenerating table visualization...")

# Column headers - include both linear and power-law fits
headers = ['L', r'$a^2$', r'$u_c$', r'$T_{50}$', 'Time Window',
           'Pearson r (P)', r'Spearman $\rho$ (P)', 
           r'Linear Slope ($R^2$)', r'Power $n$ ($R^2$)']

def format_pvalue(p):
    """Format p-value"""
    if p is None:
        return "N/A"
    if p < 0.001:
        return f"{p:.2e}"
    else:
        return f"{p:.4f}"

def format_correlation(r, p):
    """Format correlation with p-value in parentheses"""
    if r is None or p is None:
        return "N/A"
    return f"{r:.3f} ({format_pvalue(p)})"

def format_slope(slope, r2=None):
    """Format power-law exponent n with R²"""
    if slope is None:
        return "N/A"
    if r2 is not None:
        return f"{slope:.3f} ({r2:.3f})"
    return f"{slope:.3f}"

def window_name(window_type, start, end):
    """Format time window name"""
    if window_type == 'minus_2pct_full':
        return f'$T_{{50}}-0.02T_{{full}}$ ({start}-{end})'
    elif window_type == 'plus_2pct_full':
        return f'$T_{{50}}+0.02T_{{full}}$ ({start}-{end})'
    return f'{window_type} ({start}-{end})'

# Group by L
data_by_L = defaultdict(list)
for result in results:
    L_key = result['Lx'] * 10
    data_by_L[L_key].append(result)

# Generate tables for each L value
for L_value in sorted(data_by_L.keys()):
    results_L = data_by_L[L_value]
    table_data_L = []
    
    # Track correlations for potential highlighting
    pearson_values = []
    spearman_values = []
    
    # Group by (afr, threshold) to organize rows properly
    param_groups = defaultdict(list)
    for result in results_L:
        key = (result['afr'], result['threshold'])
        param_groups[key].append(result)
    
    # Sort by afr, then threshold
    for (afr, threshold) in sorted(param_groups.keys()):
        group_results = param_groups[(afr, threshold)]
        # Sort by window_type to ensure consistent order
        group_results.sort(key=lambda x: x['window_type'])
        
        for idx, result in enumerate(group_results):
            pearson_str = format_correlation(result['pearson_r'], result['pearson_p'])
            spearman_str = format_correlation(result['spearman_r'], result['spearman_p'])
            linear_str = format_slope(result['linear_slope'], result['linear_r2'])
            power_str = format_slope(result['power_exponent'], result['power_r2'])
            window_str = window_name(result['window_type'], result['step_start'], result['step_end'])
            
            # Only show L, afr, threshold, T50 in first row of each parameter group
            if idx == 0:
                table_row = [
                    str(L_value),
                    f"{afr:.2f}",
                    f"{threshold:.2f}",
                    str(result['T50']),
                    window_str,
                    pearson_str,
                    spearman_str,
                    linear_str,
                    power_str
                ]
            else:
                table_row = [
                    '',
                    '',
                    '',
                    '',
                    window_str,
                    pearson_str,
                    spearman_str,
                    linear_str,
                    power_str
                ]
            
            table_data_L.append(table_row)
            pearson_values.append(result['pearson_r'])
            spearman_values.append(result['spearman_r'])
    
    # Create figure - compact for 9 columns
    fig, ax = plt.subplots(figsize=(16, len(table_data_L) * 0.28 + 1.2))
    ax.axis('tight')
    ax.axis('off')
    
    # Add centered title addressing reviewer's question
    title_text = f'Curvature-Speed Relationship Analysis for L={L_value} μm'
    plt.title(title_text, fontsize=14, weight='bold', pad=-30, y=0.99, ha='center')
    
    # Create table
    table = ax.table(cellText=table_data_L, colLabels=headers, cellLoc='center', loc='center')
    
    # Style the table - very compact
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)
    
    # Remove all cell borders first and set compact padding
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.set_edgecolor('white')
        cell.set_facecolor('white')
        cell.PAD = 0.03
    
    # Header styling
    for i in range(len(headers)):
        cell = table[(0, i)]
        cell.set_text_props(weight='bold', color='black', fontsize=11)
        cell.set_facecolor('white')
    
    # Add three main horizontal lines (LaTeX booktabs style)
    for j in range(len(headers)):
        # Top line (above header)
        table[(0, j)].visible_edges = 'T'
        table[(0, j)].set_edgecolor('black')
        table[(0, j)].set_linewidth(1.5)
        
        # Middle line (below header)
        table[(1, j)].visible_edges = 'T'
        table[(1, j)].set_edgecolor('black')
        table[(1, j)].set_linewidth(1.0)
        
        # Bottom line (below last row)
        table[(len(table_data_L), j)].visible_edges = 'B'
        table[(len(table_data_L), j)].set_edgecolor('black')
        table[(len(table_data_L), j)].set_linewidth(1.5)
    
    # Adjust column widths for very compact layout (9 columns)
    for i in range(len(table_data_L) + 1):
        table[(i, 0)].set_width(0.05)  # L column
        table[(i, 1)].set_width(0.06)  # afr column
        table[(i, 2)].set_width(0.06)  # u_c column
        table[(i, 3)].set_width(0.07)  # T50 column
        table[(i, 4)].set_width(0.22)  # Time Window column
        table[(i, 5)].set_width(0.16)  # Pearson column
        table[(i, 6)].set_width(0.16)  # Spearman column
        table[(i, 7)].set_width(0.13)  # Linear slope column
        table[(i, 8)].set_width(0.13)  # Power exponent column
    
    # Apply red color to correlation values < 0.9
    for i in range(1, len(table_data_L) + 1):
        row_idx = i - 1
        # Pearson column (column 5)
        if pearson_values[row_idx] is not None and pearson_values[row_idx] < 0.9:
            table[(i, 5)].set_text_props(color='red')
        # Spearman column (column 6)
        if spearman_values[row_idx] is not None and spearman_values[row_idx] < 0.9:
            table[(i, 6)].set_text_props(color='red')
    
    # Add subtle lines between different parameter groups (every 2 rows)
    row_count = 0
    for (afr, threshold) in sorted(param_groups.keys()):
        n_rows = len(param_groups[(afr, threshold)])
        row_count += n_rows
        if row_count < len(table_data_L):  # Don't add line after last group
            for j in range(len(headers)):
                cell = table[(row_count + 1, j)]  # +1 because row 0 is header
                cell.visible_edges = 'T'
                cell.set_edgecolor('gray')
                cell.set_linewidth(0.5)
    
    plt.tight_layout()
    filename = f'sensitive_analysis_vk_slopes_{L_value}.png'
    plt.savefig(filename, dpi=600, bbox_inches='tight', facecolor='white', pad_inches=0.02)
    print(f"✓ Table for L={L_value} saved as '{filename}'")
    
    # Count valid correlations
    valid_pearson = sum(1 for v in pearson_values if v is not None)
    valid_spearman = sum(1 for v in spearman_values if v is not None)
    print(f"  Valid correlations: Pearson={valid_pearson}/{len(pearson_values)}, Spearman={valid_spearman}/{len(spearman_values)}")
    
    plt.close()

print("\n" + "="*80)
print("Summary:")
print(f"- Total configurations processed: {len(results)}")
print(f"- Tables generated: 3 (L=400, 600, 800)")
print(f"- Table style: LaTeX three-line (booktabs), compact")
print(f"- Font: Times New Roman")
print("="*80)
print("\nFitting Models:")
print("1. Linear fit: v = a·κ + b (direct linear relationship)")
print("   - Slope 'a' measures linear coupling strength")
print("2. Power-law fit: v = A·κ^n (via log-log regression)")
print("   - Exponent 'n' quantifies non-linear curvature-speed coupling")
print("3. Compare R² values to determine which model fits better")
print("="*80)
print("\nTime Windows:")
print("- Before T50: from T50-0.02*Tfull to T50")
print("- After T50: from T50 to T50+0.02*Tfull")
print("- Pearson & Spearman: correlation coefficients (match pair_correlations.csv)")
print("="*80)

