import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys
import re
import clr # Библиотека для работы с .NET

repeats_path = 'A_thaliana'
file_name = 'all_Ataliana_edpnew_promoters.fasta'

filter_plot_directory = 'Filter_plots'
Filter_plot_directory_path = repeats_path + "\\" + filter_plot_directory

search_patterns = []
ratios = []

# Указываем путь к папке с нашей DLL (текущая директория)
dll_path = os.path.abspath(".")
sys.path.append(dll_path)
# Загружаем библиотеку по имени файла (без расширения .dll)
clr.AddReference("DNAMotif_01")
# Импортируем класс из пространства имен C#
from DNAMotif_01 import DNAMotif, Motif


def calculate_pwm_from_text(pfm_text, background=0.25, pseudocount=0.01):
    # Парсит текстовую PFM и преобразует её в PWM (Log-Odds).
    # 1. Извлекаем строки данных (A, C, G, T)
    # Ищем символ нуклеотида и всё, что внутри скобок [ ... ]
    rows = re.findall(r'([ACGT])\s+\[\s*([\d\s]+)\]', pfm_text)
    # Создаем словарь, чтобы гарантировать порядок A, C, G, T
    pfm_dict = {char: [float(x) for x in values.split()] for char, values in rows}
    # Превращаем в 2D массив в строгом порядке ACGT
    pfm = np.array([pfm_dict[n] for n in "ACGT"])
    # 2. Добавляем псевдоотсчеты
    pfm_adj = pfm + pseudocount
    # 3. Превращаем в вероятности (PPM)
    column_sums = pfm_adj.sum(axis=0)
    ppm = pfm_adj / column_sums
    # 4. Вычисляем Log-Odds (PWM)
    pwm = np.log2(ppm / background)
    return pwm

# Данные из https://jaspar.elixir.no/matrix/MA0079.3/?revcomp=1
pfm_data = """
>TATA.1999	TBP
A [	50	835	44	892	710	848	450	358	]
C [	110	13	33	8	8	29	34	140	]
G [	45	13	9	17	5	95	164	384	]
T [	795	139	914	84	277	28	352	118	]
"""
pwm_result = calculate_pwm_from_text(pfm_data)

import itertools

octos_to_test = ['TATAAAAG',
                 'TATAAATG',
                 'TATATAAG',
                 'TTTAAAAG',
                 'TACAAAAG',

                 'TATAAGAG',
                 'AATAAAAG',
                 'TATTAAAG',
                 'TATAATAG',
                 'TATAAGAG',

                 'TATAAACG',

                 'TGTAAAAG',
                 'TCTAAAAG',
                 'TAAAAAAG',
                 'TATGAAAG',
                 'TATACAAG',]
def get_threshold_for_top_n(pwm, target_count=6144):
    """
    Находит порог (threshold) для PWM-матрицы (4xL),
    при котором ровно target_count комбинаций 8-меров проходят фильтр.
    """
    nucs = ['A', 'C', 'G', 'T']
    motif_len = len(pwm) // 4
    process_motif = DNAMotif()
    process_motif.SetMatrix(pwm)
    process_motif.SetThreshold(-1e6)
    all_scores = []

    # Генерируем все 65536 комбинаций (для 8-ми нуклеотидов)
    for combo in itertools.product(nucs, repeat=motif_len):
        seq_str = "".join(combo)
        # Вызываем ваш стандартный вычислитель
        found_motifs = process_motif.FindMotifsTATA1999(seq_str)
        # Берем вес найденного мотива (он там будет один, так как длина совпадает)
        if found_motifs:
            all_scores.append(found_motifs[0].Weight)

    # Сортируем веса от самого высокого к самому низкому
    print(f"Array length: {len(all_scores)}")
    all_scores.sort(reverse=True)
    # Берем значение веса на позиции target_count
    # Если их несколько с одинаковым весом, все они пройдут при этом пороге
    if target_count <= len(all_scores):
        calculated_threshold = all_scores[target_count - 1]
    else:
        calculated_threshold = all_scores[-1]

    for oct in octos_to_test:
        found_motifs = process_motif.FindMotifsTATA1999(oct)
        # Берем вес найденного мотива (он там будет один, так как длина совпадает)
        if found_motifs:
            print(f"Weight for {oct} - {found_motifs[0].Weight:.4f}")


    print(f"Calculated threshold for top 6144: {calculated_threshold:.4f}")
    return calculated_threshold


def process_fasta_with_pwm(pwm, threshold=0.0):
    global repeats_path, file_name
    file_path = os.path.join(repeats_path, file_name)

    total_sequences = 0
    no_pattern_count = 0
    motif_len = len(pwm) // 4
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Списки для разделения последовательностей
    TATA = []
    TATAless = []

    process_motif = DNAMotif()
    process_motif.SetMatrix(pwm)
    process_motif.SetThreshold(threshold)

    # Инициализация счетчика позиций (берем длину из первой записи)
    with open(file_path, "r") as handle:
        records = SeqIO.parse(handle, "fasta")
        for record in records:
            req_len = len(str(record.seq)) - motif_len + 1
            position_counts = [0] * req_len
            break

        for record in records:
            total_sequences += 1
            seq_str = str(record.seq).upper()  # Convert to uppercase
            found_in_record = False
            max_score = -1e10

            found_motifs = process_motif.FindMotifsTATA1999(seq_str)
            for motif in found_motifs:
                if 0 <= motif.Position < len(position_counts):
                    # position_counts[motif.Position] += 0.5 * (motif.Weight - threshold)
                    position_counts[motif.Position] += 1
                    found_in_record = True
                if motif.Weight > max_score:
                    max_score = motif.Weight

            if (len(found_motifs) > 0):
                TATA.append(record)
            else:
                TATAless.append(record)
                no_pattern_count += 1
            # print (np.round(max_score, 2), end=', ')

    file_name_list = file_name.split('_')
    file_name_list[0] = 'TATA2'
    new_file_name = '_'.join(file_name_list)
    output_fasta = os.path.join(repeats_path, new_file_name)
    with open(output_fasta, "w") as out_handle:
        for record in TATA:
            SeqIO.write([record], out_handle, "fasta")

    file_name_list[0] = 'TATA2less'
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
    motif_len = len(pwm_lists) // 4
    init = -50
    end = 30
    x_positions = list(range(init, end - motif_len + 1))  # [-50, -49, ..., 0 (TSS), ..., +35]

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
    plt.margins(x=0)
    plt.tight_layout()

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
    pwm_matrix = calculate_pwm_from_text(pfm_data)
    pwm_lists = pwm_matrix.flatten().tolist()
    # ВЫЧИСЛЯЕМ ПОРОГ:
    # threshold = get_threshold_for_top_n(pwm_lists, 6144)
    # exit(0)
    threshold = 0.7
    # print(*search_patterns, sep="\n")
    full_filepath = os.path.join(repeats_path, file_name)
    matrix_name = pfm_data.split('\t')[0][2:]
    matrix_name = matrix_name.replace('.', '_')
    matrix_str = matrix_name + f" threshold = {threshold:.2f}"

    fn_list = file_name.split('_')
    output_file = fn_list[0] + '_' + fn_list[1] + "_TATAmat" + matrix_str + " .txt"
    output_file = os.path.join(Filter_plot_directory_path, output_file)
    # Main processing logic
    frequencies, total_sequences, no_pattern_count = process_fasta_with_pwm(pwm_lists, threshold)
    save_results_to_file(output_file, frequencies, total_sequences, no_pattern_count)

    output_file = fn_list[0] + '_' + fn_list[1] + '_TATAmat' + matrix_str + '.png'
    output_file = os.path.join(Filter_plot_directory_path, output_file)
    prefix = '. Percent of TATAmat boxes. ' + matrix_str
    plot_nucleotides(ratios, output_file, prefix)
