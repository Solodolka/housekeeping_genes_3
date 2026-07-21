import os
from Bio import SeqIO
import matplotlib.pyplot as plt

# Possible variants of the TATAWAWR sequence
pattern = "TATAWAWR"
# pattern = "SSRCGCC"
# pattern = "RTDKKKK"
pattern_name = "TATAWAWR"
nuc_num = len(pattern)  # Length of the TATAWAWR sequence

repeats_path = os.getenv('MY_GLOBAL_VAR')
if repeats_path is None:
    repeats_path = 'H_vulgare'

organizm = repeats_path.replace('_', '')
file_name = f'all_{organizm}_epdnew_promoters.fasta'

filter_plot_directory = 'Filter_plots'
Filter_plot_directory_path = repeats_path + "\\" + filter_plot_directory

search_patterns = []
ratios = []

def expand_pattern():
    global pattern, search_patterns
    # Словарь замен
    iupac = {
        'A': 'A',
        'G': 'G',
        'C': 'C',
        'T': 'T',
        'W': 'AT',
        'R': 'AG',
        'S': 'GC',
        'D': 'AGT',
        'K': 'GT',
    }
    search_patterns = []

    def solve(current_seq, index):
        # Если дошли до конца паттерна, сохраняем готовую строку
        if index == len(pattern):
            search_patterns.append(current_seq)
            return
        # Берем текущий символ и смотрим, какие у него есть варианты
        char = pattern[index]
        options = iupac[char]
        # Перебираем все варианты для данной позиции
        for nucleotide in options:
            solve(current_seq + nucleotide, index + 1)

    solve("", 0)
    return search_patterns


def find_positions_in_record(seq):
    """
    Return a list of positions where any of the searched patterns are found.
    :param seq: DNA sequence string
    :return: List of integers representing starting indices of patterns
    """
    positions = []
    for pattern in search_patterns:
        pos = seq.find(pattern)
        while pos != -1:
            positions.append(pos)
            pos = seq.find(pattern, pos + 1)
    return positions


def process_fasta_file():
    global repeats_path, file_name
    """
    Process the entire FASTA file and return statistics on pattern distribution.
    :param file_path: Full path to the FASTA file
    :return: Tuple (dictionary of positions and their frequency, total number of sequences, number of promoters without pattern)
    """
    total_sequences = 0
    no_pattern_count = 0
    file_path = os.path.join(repeats_path, file_name)

    TATA = []
    TATAless = []
    with open(file_path, "r") as handle:
        records = SeqIO.parse(handle, "fasta")
        for record in records:
            position_counts = [0] * (len(str(record.seq))-nuc_num+1)
            break

        for record in records:
            total_sequences += 1
            seq_str = str(record.seq).upper()  # Convert to uppercase
            # Check presence of any patterns
            found_positions = find_positions_in_record(seq_str)
            if not found_positions:
                no_pattern_count += 1
            else:
                for pos in found_positions:
                    position_counts[pos] += 1

            if (found_positions):
                TATA.append(record)
            else:
                TATAless.append(record)


    file_name_list = file_name.split('_')
    file_name_list[0] = pattern_name
    new_file_name = '_'.join(file_name_list)
    output_fasta = os.path.join(repeats_path, new_file_name)
    with open(output_fasta, "w") as out_handle:
        for record in TATA:
            SeqIO.write([record], out_handle, "fasta")

    file_name_list[0] = f'{pattern_name}_less'
    new_file_name = '_'.join(file_name_list)
    output_fasta = os.path.join(repeats_path, new_file_name)
    with open(output_fasta, "w") as out_handle:
        for record in TATAless:
            SeqIO.write([record], out_handle, "fasta")

    return position_counts, total_sequences, no_pattern_count


def save_results_to_file(output_filename, position_freqs, total_seq_count, no_pattern_count):
    global ratios
    """
    Save analysis results into a text file.
    :param output_filename: Name of the output file
    :param position_freqs: Dictionary mapping positions to their frequencies
    :param total_seq_count: Total number of processed sequences
    :param no_pattern_count: Number of promoters without the target pattern
    """
    with open(output_filename, 'w', encoding='utf-8') as f:
        # General information header
        f.write(f'# Total sequences: {total_seq_count}\n')
        f.write(f'# Promoters without pattern: {no_pattern_count}\n\n')
        f.write(f'# Share without pattern: {no_pattern_count/total_seq_count}\n\n')

        # Table with detailed position-based data
        f.write('Position\tFrequency\tRelative share (%)\n')
        ratios = []
        for i, freq in enumerate(position_freqs):  # Include all possible positions up to maximum one found
            ratio = freq / total_seq_count * 100 if total_seq_count > 0 else 0
            ratios.append(ratio)
            f.write(f'{i}\t{freq}\t{ratio:.2f}%\n')

        print("File successfully written:", output_filename)


def plot_nucleotides(data, output_filename, prefix):
    # Создание массива осей X с новыми индексами
    x_positions = list(range(-50, 23))  # [-50, -49, ..., 0 (TSS), ..., +35]

    # Построение графика
    plt.figure(figsize=(10, 3))
    plt.plot(x_positions, data[:len(x_positions)], marker='.', markersize=0, color='#00008f')

    # Установка собственных меток для оси X
    step = 10  # Метки ставятся через каждые 10 единиц
    custom_ticks = x_positions[::step]  # Берём каждую десятую позицию
    custom_labels = [str(i) if i != 0 else 'TSS' for i in custom_ticks]
    plt.xticks(custom_ticks, custom_labels)  # Применяем новые метки и позиции

    # Подписываем оси
    plt.axvline(x=0, color="red", linestyle="--", linewidth=0.5)  # Линия вертикальной разметки
    plt.axvline(x=-28, color="blue", linestyle="--", linewidth=0.5)  # Линия вертикальной разметки

    fn_list = file_name.split('_')
    str0 = fn_list[0][0].upper() + fn_list[0][1:]
    str1 = f"{fn_list[1][0].capitalize()}.{fn_list[1][1:]}"

    # Настройка внешнего вида графика
    title = str0 + ' ' + str1
    title += prefix
    plt.title(title)
    # plt.ylabel()
    plt.grid(True)

    # Сохраняем график
    plt.savefig(output_filename)
    plt.close()
    print(f'Graph saved as {output_filename}')


def create_directory_if_not_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f'Directory "{dir_path}" created.')
    else:
        print(f'Directory "{dir_path}" already exists.')


if __name__ == "__main__":
    create_directory_if_not_exists(Filter_plot_directory_path)
    # Define paths
    expand_pattern()
    # print(*search_patterns, sep="\n")
    full_filepath = os.path.join(repeats_path, file_name)

    fn_list = file_name.split('_')
    output_file = fn_list[0] + '_' + fn_list[1] + f"_{pattern_name}_box.txt"
    output_file = os.path.join(Filter_plot_directory_path, output_file)
    # Main processing logic
    frequencies, total_sequences, no_pattern_count = process_fasta_file()
    save_results_to_file(output_file, frequencies, total_sequences, no_pattern_count)

    output_file = fn_list[0] + '_' + fn_list[1] + f"_{pattern_name}_box.png"
    output_file = os.path.join(Filter_plot_directory_path, output_file)
    prefix = f'. Percent of {pattern} boxes'
    plot_nucleotides(ratios, output_file, prefix)
