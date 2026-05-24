import os
from Bio import SeqIO

# Константы
repeats_path = 'HK_promoters'
file_name = 'unfiltered_barley_epdnew_hv_promoters_all.fasta'
results_path = 'Nucleotides_probabilities'
SHIFTS_NUM = 90  # Число позиций (нуклеотиды), которые будут обработаны
RUN_NUMBER = 1  # Текущий номер прохода файла
COUNT_BEG_VALUES = [[0] * 4 for _ in range(SHIFTS_NUM)]  # Массив накопленных сумм A,C,T,G в позициях начала
Nucleotides = ['A', 'C', 'T', 'G']  # Возможные типы нуклеотидов

def read_and_process_repeats(directory):
    """Перебираем все файлы с повторами"""
    global RUN_NUMBER
    filepath = os.path.join(directory, file_name)
    for record in SeqIO.parse(filepath, "fasta"):
        process_repeat_beg(record.seq)  # Обрабатываем начало
        RUN_NUMBER += 1

def process_repeat_beg(repeat):
    """Обрабатываем начало последовательности"""
    values = [[0] * 4 for _ in range(SHIFTS_NUM)]

    for pos in range(SHIFTS_NUM):
        nucleotide = repeat[pos]
        if nucleotide in Nucleotides:
            idx = Nucleotides.index(nucleotide)
            values[pos][idx] = 1  # Присваиваем единицу соответствующему типу нуклеотида

    add_to_count(values, COUNT_BEG_VALUES)

def add_to_count(values, count_matrix):
    """Накопление результатов"""
    for row_idx, row in enumerate(values):
        for col_idx, val in enumerate(row):
            count_matrix[row_idx][col_idx] += val

def normalize_counts(count_matrix):
    """Нормализация полученных сумм в вероятности"""
    norm_matrix = []
    for row in count_matrix:
        total = sum(row)
        prob_row = [val / total for val in row]
        norm_matrix.append(prob_row)
    return norm_matrix

def write_results(values):
    """Запись результата в файл"""
    prefix = file_name.split('_')[1]
    filename = os.path.join(results_path, 'result_' + prefix + '.txt')
    with open(filename, 'w', encoding='utf-8') as f:
        for row in values:
            formatted_row = "\t".join([f"{val:.8f}" for val in row]).replace(".", ".")
            f.write(formatted_row + "\n")

read_and_process_repeats(repeats_path)
normalized_beg_probs = normalize_counts(COUNT_BEG_VALUES)
write_results(normalized_beg_probs)
