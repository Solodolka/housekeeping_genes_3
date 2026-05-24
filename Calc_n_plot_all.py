from Bio import SeqIO
import matplotlib.pyplot as plt
from weblogo import *

repeats_path = 'A_thaliana'
organizm = repeats_path.replace('_', '')

file_name_list = [
    f'all_{organizm}_epdnew_promoters.fasta',
    f'NoProhibNuc_{organizm}_epdnew_promoters.fasta',
    f'ProhibNuc_{organizm}_epdnew_promoters.fasta',
    # f'TATA2_{organizm}_epdnew_promoters.fasta',
    # f'TATA2less_{organizm}_epdnew_promoters.fasta',
]

prob_plot_directory = 'Probabilities_plots'
Prob_plot_directory_path = repeats_path + "\\" + prob_plot_directory
NUCLEOTIDES_IN_STRING = 80  # Number of positions (nucleotides) to process
COUNT_BEG_VALUES = [[0] * 4 for _ in range(NUCLEOTIDES_IN_STRING)]  # Matrix for accumulated sums of A,C,T,G in start positions
Nucleotides = ['A', 'C', 'T', 'G']  # Possible nucleotide types


def read_and_process_repeats(file_name):
    """Iterating through all sequences from FASTA file."""
    filepath = os.path.join(repeats_path, file_name)
    for record in SeqIO.parse(filepath, "fasta"):
        if len(record) < NUCLEOTIDES_IN_STRING:
            print(f"Длина записи другая! {record}. Длина {len(record)} ")
            exit(1)
        process_repeat_beg(record.seq)  # Process beginning of sequence

def process_repeat_beg(repeat):
    """Processing start of sequence."""
    values = [[0] * 4 for _ in range(NUCLEOTIDES_IN_STRING)]

    for pos in range(NUCLEOTIDES_IN_STRING):
        nucleotide = repeat[pos].upper()  # Convert to uppercase for reliability
        if nucleotide in Nucleotides:
            idx = Nucleotides.index(nucleotide)
            values[pos][idx] = 1  # Assign value 1 to corresponding nucleotide type

    add_to_count(values, COUNT_BEG_VALUES)

def add_to_count(values, count_matrix):
    """Accumulate counts."""
    for row_idx, row in enumerate(values):
        for col_idx, val in enumerate(row):
            count_matrix[row_idx][col_idx] += val

def normalize_counts(count_matrix):
    """Normalize sums to probability values."""
    norm_matrix = []
    for row in count_matrix:
        total = sum(row)
        prob_row = [val / total for val in row]
        norm_matrix.append(prob_row)
    return norm_matrix


def plot_nucleotides(norm_values, file_name):
    """
    Plot graphs representing probabilities of four nucleotide types using normalized data.
    """
    positions = list(range(-50, NUCLEOTIDES_IN_STRING - 50))
    probabilities = {'A': [], 'C': [], 'G': [], 'T': []}

    # Distribute normalized values among respective nucleotides
    for row in norm_values:
        probabilities['A'].append(row[0])
        probabilities['C'].append(row[1])
        probabilities['T'].append(row[2])
        probabilities['G'].append(row[3])

    # Plot graphs
    plt.figure(figsize=(10, 6))
    colors = ['black', 'blue', 'yellow', 'red']
    labels = ['Adenine (A)', 'Cytosine (C)', 'Guanine (G)', 'Thymine (T)']

    for nucleotide, color, label in zip(probabilities.keys(), colors, labels):
        plt.plot(positions, probabilities[nucleotide], marker='.', markersize=0, color=color, label=label)

    # Custom x-axis ticks setup
    step = 10  # Step between marks on x-axis
    custom_ticks = positions[::step]
    custom_labels = [str(i) if i != 0 else 'TSS' for i in custom_ticks]
    plt.xticks(custom_ticks, custom_labels)

    # Additional vertical lines
    plt.axvline(x=0, color="red", linestyle="--", linewidth=0.5)
    plt.axvline(x=-28, color="blue", linestyle="--", linewidth=0.5)

    fn_list = file_name.split('_')
    str0 = fn_list[0][0].upper() + fn_list[0][1:]
    str1 = f"{fn_list[1][0].capitalize()}.{fn_list[1][1:]}"

    # Graph layout
    plt.legend(loc="upper right")
    plt.title(str0 + ' ' + str1 +'. Probability of nucleotide presence by position')
    plt.xlabel('Position')
    plt.ylabel('Probability')
    plt.grid(True)

    # Save graph image
    output_filename = fn_list[0] + '_' + fn_list[1] + '_prob_graphs.png'
    output_filename = os.path.join(Prob_plot_directory_path, output_filename)
    plt.savefig(output_filename)
    plt.close()
    print(f'Graph saved as {output_filename}')

def create_weblogo(file_name):
    filepath = os.path.join(repeats_path, file_name)
    fin = open(filepath, "r")
    seqs = read_seq_data(fin)

    fn_list = file_name.split('_')
    str0 = fn_list[0][0].upper() + fn_list[0][1:]
    str1 = f"{fn_list[1][0].capitalize()}.{fn_list[1][1:]}"

    # 2. Создание логотипа
    logo = LogoData.from_seqs(seqs)
    options = LogoOptions()
    options.logo_title = str0 + ' ' + str1 + " WebLogo"
    # Настройка параметров логотипа
    # делаем столько стопок, сколько нужно, чтобы занять одну строку
    options.stacks_per_line = NUCLEOTIDES_IN_STRING + 1
    options.stack_width = std_sizes['medium']
    options.stack_aspect_ratio = 30
    options.color_scheme = ColorScheme([
        SymbolColor("G", "yellow"),  # Guanine (Гуанин) зеленым
        SymbolColor("C", "blue"),  # Cytosine (Цитозин) желтым
        SymbolColor("A", "black"),  # Adenine (Аденин) красным
        SymbolColor("T", "red")  # Thymine (Тимин) синим
    ])
    options.yaxis_scale = 1.0
    options.yaxis_tic_interval = 0.1
    options.fineprint_font_size = 30
    options.logo_height = 5
    options.resolution = 600
    options.unit_name = "bits"
    options.stack_spacing = 0  # Минимизирует расстояние между стопками
    options.show_fineprint = True  # Отключает подпись снизу
    options.logo_margin = 30
    format = LogoFormat(logo, options)

    # 3. Save graph image
    image_bytes = jpeg_formatter(logo, format)
    output_filename = fn_list[0] + '_' + fn_list[1] + '_weblogo.jpg'
    output_filename = os.path.join(Prob_plot_directory_path, output_filename)
    with open(output_filename, "wb") as f:
        f.write(image_bytes)
    print(f'Graph saved as {output_filename}')

def create_directory_if_not_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f'Directory "{dir_path}" created.')
    else:
        print(f'Directory "{dir_path}" already exists.')

if __name__ == "__main__":
    for file_name in file_name_list:
        create_directory_if_not_exists(Prob_plot_directory_path)
        # Processing sequences
        read_and_process_repeats(file_name)
        # Normalize collected sums to probabilities
        normalized_probs = normalize_counts(COUNT_BEG_VALUES)
        # Immediately build and save the plot
        plot_nucleotides(normalized_probs, file_name)
        # Create WebLogo
        create_weblogo(file_name)