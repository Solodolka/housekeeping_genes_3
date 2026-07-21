import os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
import sys
import re
import clr # Библиотека для работы с .NET


repeats_path = os.getenv('MY_GLOBAL_VAR')
if repeats_path is None:
    repeats_path = 'M_musculus'

organizm = repeats_path.replace('_', '')
file_name = f'TATAless_{organizm}_epdnew_promoters.fasta'
ratios = []

# Указываем путь к папке с нашей DLL (текущая директория)
dll_path = os.path.abspath(".")
sys.path.append(dll_path)
# Загружаем библиотеку по имени файла (без расширения .dll)
clr.AddReference("DNAMotif_01")
# Импортируем класс из пространства имен C#
from DNAMotif_01 import DNAMotif, Motif

def calculate_pwm_from_text(pfm_text, background=0.25, pseudocount=0.1):
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
>MA0079.2	SP1
A  [     0      0      0      4      2      0      1      0      6      3 ]
C  [    32     30     35     27      5     28     31     24     25     26 ]
G  [     1      1      0      0     15      1      0      3      0      3 ]
T  [     2      4      0      4     13      6      3      8      4      3 ]
"""

pwm_result = calculate_pwm_from_text(pfm_data)

def score_sequence(window, pwm, nuc_to_idx):
    """Вычисляет score для конкретного окна последовательности."""
    score = 0
    for i, nuc in enumerate(window):
        if nuc in nuc_to_idx:
            score += pwm[nuc_to_idx[nuc], i]
        else:
            score += 0  # Для N или неизвестных символов
    return score


def process_fasta_with_pwm(pwm, threshold=0.0):
    global repeats_path, file_name
    file_path = os.path.join(repeats_path, file_name)

    total_sequences = 0
    no_pattern_count = 0
    motif_len = len(pwm) // 4
    nuc_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    # Списки для разделения последовательностей
    FoundMotif = []
    NoMotif = []

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

            found_motifs = process_motif.FindMotifs(seq_str)
            for motif in found_motifs:
                if 0 <= motif.Position < len(position_counts):
                    position_counts[motif.Position] += 0.5 * (motif.Weight - threshold)
                    found_in_record = True
                if motif.Weight > max_score:
                    max_score = motif.Weight

            if found_in_record:
                ...
                # FoundMotif.append(record)
            else:
                max_score = 0
                no_pattern_count += 1
                # NoMotif.append(record)
            # print (np.round(max_score, 2), end=', ')

    return position_counts, total_sequences, no_pattern_count


def plot_nucleotides(data, output_filename, prefix):
    # 1. Определяем длину данных
    data_len = len(data)
    if data_len == 0:
        print("No data to plot.")
        return

    # 2. Создаем массив X, соответствующий длине данных.
    #    Позиция 0 в массиве data соответствует -28 п.н. от TSS.
    #    Позиция data_len-1 соответствует +35 п.н. от TSS.
    #    Массив X будет от -28 до -28 + data_len - 1
    x_positions = list(range(-50, -50 + data_len))

    # Построение графика
    plt.figure(figsize=(15, 4.5))
    plt.plot(x_positions, data, marker='.', markersize=0, color='#00008f')

    # 3. Настройка меток на оси X
    # Шаг меток: ставим их через каждые 5-10 позиций, чтобы они не сливались.
    step = 5
    custom_ticks = x_positions[::step]
    custom_labels = [str(i) if i != 0 else 'TSS' for i in custom_ticks]
    plt.xticks(custom_ticks, custom_labels)

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


if __name__ == "__main__":
    pwm_matrix = calculate_pwm_from_text(pfm_data)
    pwm_lists = pwm_matrix.flatten().tolist()

    full_filepath = os.path.join(repeats_path, file_name)

    fn_list = file_name.split('_')
    output_file = fn_list[0] + '_' + fn_list[1] + f"_GC_box.txt"
    output_file = os.path.join(repeats_path, output_file)
    # Main processing logic
    frequencies, total_sequences, no_pattern_count = process_fasta_with_pwm(pwm_lists)
    save_results_to_file(output_file, frequencies, total_sequences, no_pattern_count)

    output_file = fn_list[0] + '_' + fn_list[1] + f"_GC_box.png"
    output_file = os.path.join(repeats_path, output_file)
    prefix = '. Percent of GC boxes'
    plot_nucleotides(ratios, output_file, prefix)
