import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import os

# Create output directory if it doesn't exist
os.makedirs('Publication_figures', exist_ok=True)

# Data embedded directly in the script
group1_data = """Bacillus_subtilis	1-1-1-1-3-2-1
Enterococcus_faecalis	2-1-1
Escherichia_coli	2-2-1-1-1
Lactobacillus_fermentum	1-2-1-1
Listeria_monocytogenes	1-2-1-1-1
Pseudomonas_aeruginosa	4
Salmonella_enterica	3-1-2-1
Staphylococcus_aureus	1-1-1-1-1-1"""

group2_data = """Acinetobacter_baumannii	6
Bacillus_pacificus	1
Bacteroides_vulgatus	1
Bifidobacterium_adolescentis	3-1-1
Clostridium_beijerinckii	1-1-1-1-1-1-1-1-1-1-1-1-1-1
Cutibacterium_acnes	2
Deinococcus_radiodurans	1-2
Enterococcus_faecalis	3-1
Escherichia_coli	1-1-1-1-2-1
Helicobacter_pylori	1-1
Lactobacillus_gasseri	1
Neisseria_meningitidis	4
Porphyromonas_gingivalis	4
Pseudomonas_aeruginosa	1
Rhodobacter_sphaeroides	4
Schaalia_odontolytica	2
Staphylococcus_aureus	1-2-1-1-1
Staphylococcus_epidermidis	1-1-1-2-1
Streptococcus_agalactiae	7
Streptococcus_mutans	2-3"""

group3_data = """Acetivibrio_thermocellus	2-1-1
Allomeiothermus_silvanus	2
Clostridium_perfringens	1-1-1-1-1-1-1-1
Coraliomargarita_akajimensis	2
Corynebacterium_glutamicum	1-1-1-2-1
Desulfoscipio_gibsoniae	1-1-1-1-1-1-1-1
Desulfosporosinus_acidiphilus	1-1-1-1-1-1-1-1-1
Desulfosporosinus_meridiei	1-1-2-1-1-1-1-1-1-1
Echinicola_vietnamensis	4
Escherichia_coli	1-1-1-1-1-2
Fervidobacterium_pennivorans	1
Frateuria_aurantia	4
Hirschia_baltica	1-1
Nocardiopsis_dassonvillei	2-1-1-1
Olsenella_uli	1
Salmonella_bongori	1-2-2-1-1
Salmonella_enterica	2-1-1-1-1-1
Sediminispirochaeta_smaragdinae	1-1
Segniliparus_rotundus	1
Streptococcus_pyogenes	2-3-1
Stutzerimonas_stutzeri	1-3
Terriglobus_roseus	2
Thermobacillus_composti	2-1-1-1"""

def format_species_name(species):
    """Convert species name to italicized format with first initial and full species name."""
    parts = species.split('_')
    if len(parts) >= 2:
        genus = parts[0]
        species_name = parts[1]
        return f"${genus[0]}. {species_name}$"  # LaTeX formatting for italics
    return species

def parse_pattern(pattern):
    """Parse the pattern string into a list of integers, sorted descending."""
    counts = [int(x) for x in pattern.split('-')]
    counts.sort(reverse=True)  # Sort descending so largest counts appear first
    return counts

def get_pattern_string(pattern_list):
    """Convert a parsed pattern list to a standardized string for sorting."""
    return '-'.join(str(x) for x in pattern_list)

def parse_data_string(data_string):
    """Parse the multi-line data string into a list of (species, pattern) tuples."""
    species_data = []
    for line in data_string.strip().split('\n'):
        parts = line.strip().split('\t')
        if len(parts) == 2:
            species, pattern = parts
            species_data.append((species, parse_pattern(pattern)))
    return species_data

# Parse all three groups
all_data = [
    parse_data_string(group1_data),
    parse_data_string(group2_data),
    parse_data_string(group3_data)
]

# Color scheme for gradient shading
def get_color_for_count(count, max_count=7):
    """Get color based on count value - darker for higher counts."""
    # Using a blue gradient from light to dark - slightly darker overall
    colors = {
        1: '#A1C6E7',  # Light blue (darker than before)
        2: '#7AB3E0',  # Light-medium blue
        3: '#529FD8',  # Medium blue
        4: '#2B8CCF',  # Medium-dark blue
        5: '#1A76C0',  # Dark blue
        6: '#0960B0',  # Darker blue
        7: '#004A9F'   # Very dark blue (darker than before)
    }
    return colors.get(count, colors[7])  # Default to darkest for any count > 7

panel_titles = ['Zymo/Titan-1', 'ATCC/16S', 'Phylotag/16S']
file_names = ['zymo_d6300_titan', 'atcc_msa_1003_16s', 'phylotag_16s']

# Maximum number of species to show per panel
max_species_per_panel = 25

# Create separate figures
for idx, (species_data, title, filename) in enumerate(zip(all_data, panel_titles, file_names)):
    # Enhanced sorting: first by total ASV count (descending), 
    # then by standardized pattern string (lexicographically)
    species_data.sort(key=lambda x: (-sum(x[1]), get_pattern_string(x[1])))
    
    # Limit the number of species shown if needed
    species_data = species_data[:max_species_per_panel]
    
    # Create a new figure for each dataset with adjusted size for better proportions
    fig, ax = plt.subplots(1, 1, figsize=(6, len(species_data) * 0.4 + 1.5))
    
    # Set the title - single line format with better positioning
    fig.suptitle(f'Amplitypes: {title}', fontsize=16, fontweight='bold', 
                 x=0.55, y=0.97)  # Centered over the bars with more space
    
    # Set up the plot
    ax.set_xlim(-3.5, 9)
    ax.set_ylim(-0.5, len(species_data) - 0.5)  # Tighter vertical spacing
    
    # Remove axes
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Add visual grouping lines for species with same total ASV count
    current_total = None
    for i, (species, pattern) in enumerate(species_data):
        y_pos = len(species_data) - i - 1
        total_asvs = sum(pattern)
        
        # Add a subtle horizontal line when total ASV count changes
        if current_total is not None and total_asvs != current_total:
            ax.axhline(y=y_pos + 0.5, color='lightgray', linestyle='-', linewidth=0.5, alpha=0.5)
        current_total = total_asvs
        
        # Add species name (italicized)
        formatted_name = format_species_name(species)
        ax.text(-0.2, y_pos, formatted_name, 
               ha='right', va='center', fontsize=9, style='italic')
        
        # Calculate positions for ASV blocks
        total_width = 8.0  # Total width for pattern display
        
        # Create blocks for each ASV count
        x_pos = 0.5
        for j, count in enumerate(pattern):
            # Calculate width proportional to ASV count
            width = (count / total_asvs) * total_width if total_asvs > 0 else 0
            
            # Create blocks for each ASV count
            if width > 0:
                # Get color based on count value
                block_color = get_color_for_count(count)
                
                # Create fancy box with rounded corners - gradient colors
                box = FancyBboxPatch((x_pos, y_pos - 0.35), width - 0.08, 0.7,
                                   boxstyle="round,pad=0.02",
                                   facecolor=block_color,
                                   edgecolor='white',
                                   linewidth=1.5)
                ax.add_patch(box)
                
                # Add number in the center of the box
                ax.text(x_pos + (width - 0.08)/2, y_pos, str(count),
                       ha='center', va='center', fontsize=9,
                       color='white', fontweight='bold')
                
                x_pos += width  # Move to next position

    # Save the figure
    plt.tight_layout()
    plt.savefig(f'Publication_figures/{filename}_amplitypes.png', dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()  # Close the figure to free memory


