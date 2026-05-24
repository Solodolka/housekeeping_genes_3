import os
from collections import Counter, defaultdict
from Bio import SeqIO

# Пути и имена файлов
repeats_path = 'S_pombe'
file_name = 'all_Spombe_epdnew_promoters.fasta'
results_path = repeats_path + '\\' + 'Probabilities_plots'

# Диапазон позиций (20-27 включительно)
START_POS = 19
END_POS = 26
POS_RANGE = END_POS - START_POS + 1  # Должно быть 8

def calc_friquancies():
    global output_file
    # Счётчики для каждой позиции и нуклеотида
    position_counters = [Counter() for _ in range(POS_RANGE)]
    total_sequences = 0

    # Читаем FASTA-файл
    fasta_file = os.path.join(repeats_path, file_name)
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).upper()
        # Проверяем, что последовательность достаточно длинная
        if len(seq) <= END_POS:
            continue
        total_sequences += 1
        # Для каждой позиции в диапазоне 20-27
        for i in range(POS_RANGE):
            pos = START_POS + i
            nucleotide = seq[pos]
            position_counters[i][nucleotide] += 1

    # Записываем результаты в файл
    fn_list = file_name.split('_')
    output_filename = fn_list[0] + '_' + fn_list[1] + '_nucleotide_probabilities.txt'
    output_filename = os.path.join(results_path, output_filename)
    with open(output_filename, 'w', encoding='utf-8') as out:
        # Заголовок
        out.write('Позиция\tA (%)\tC (%)\tG (%)\tT (%)\n')
        # Для каждой позиции
        for i in range(POS_RANGE):
            pos = START_POS + i - 50
            counter = position_counters[i]
            total = sum(counter.values())
            # Если не было последовательностей — пропускаем
            if total == 0:
                continue
            # Вычисляем проценты
            a_pct = counter['A'] / total * 100 if 'A' in counter else 0
            c_pct = counter['C'] / total * 100 if 'C' in counter else 0
            g_pct = counter['G'] / total * 100 if 'G' in counter else 0
            t_pct = counter['T'] / total * 100 if 'T' in counter else 0
            # Записываем строку
            out.write(f'{pos}\t{a_pct:.2f}\t{c_pct:.2f}\t{g_pct:.2f}\t{t_pct:.2f}\n')

calc_friquancies()