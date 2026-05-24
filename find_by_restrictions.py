import os
from Bio import SeqIO

repeats_path = 'H_vulgare'
organizm = repeats_path.replace('_', '')
file_name = f'all_{organizm}_epdnew_promoters.fasta'
filter_plot_directory = 'Filter_plots'
Filter_plot_directory_path = os.path.join(repeats_path, filter_plot_directory)

ratios = []

def check_restricted_nuc(seq):
    global repeats_path
    """
    Проверяет, соответствует ли фрагмент с 21 по 28 позицию (включительно)
    следующим условиям:
    - 2-я позиция (22-я в общей) — G или C.
    - 3-я позиция (23-я в общей) — G.
    - 4-я или 5-я позиция (24-я или 25-я в общей) — G или C.
    """

    # Получаем фрагмент с 21 по 28 позицию (индексы с 20 по 27)
    fragment = seq[20:28].upper() # for H_sapiens
    if (repeats_path == "A_thaliana"):
        fragment = seq[16:24].upper()
    if (repeats_path == "H_vulgare"):
        fragment = seq[16:24].upper()
    if (repeats_path == "Z_mays"):
        fragment = seq[16:24].upper()
    if (repeats_path == "M_mulatta"):
        fragment = seq[19:27].upper()
    if (repeats_path == "M_musculus"):
        fragment = seq[20:28].upper()

    # Проверяем условия
    if fragment[1] in ('G', 'C'):
        return True
    if fragment[2] == 'G':
        return True
    if fragment[3] in ('G', 'C'):
        return True
    if fragment[4] in ('G', 'C'):
        return True

    return False

def process_fasta_file():
    global repeats_path, file_name
    """
    Обрабатывает FASTA-файл и возвращает статистику по распределению позиций,
    а также разделяет последовательности на соответствующие и не соответствующие критерию.
    """
    total_sequences = 0
    no_pattern_count = 0
    file_path = os.path.join(repeats_path, file_name)

    matching_records = []
    non_matching_records = []

    with open(file_path, "r") as handle:
        records = list(SeqIO.parse(handle, "fasta"))  # Преобразуем в список для повторного прохода

        for record in records:
            total_sequences += 1
            seq_str = str(record.seq).upper()

            if check_restricted_nuc(seq_str):
                matching_records.append(record)
            else:
                non_matching_records.append(record)
                no_pattern_count += 1

    # Сохраняем результаты в новые FASTA-файлы
    file_name_list = file_name.split('_')
    file_name_list[0] = 'ProhibNuc'
    new_file_name = '_'.join(file_name_list)
    output_fasta_match = os.path.join(repeats_path, new_file_name)

    file_name_list[0] = 'NoProhibNuc'
    new_file_name = '_'.join(file_name_list)
    output_fasta_nomatch = os.path.join(repeats_path, new_file_name)

    with open(output_fasta_match, "w") as out_handle:
        SeqIO.write(matching_records, out_handle, "fasta")

    with open(output_fasta_nomatch, "w") as out_handle:
        SeqIO.write(non_matching_records, out_handle, "fasta")

    return total_sequences, no_pattern_count


def save_results_to_file(output_filename, total_seq_count, no_pattern_count):
    global ratios
    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write(f'# Total sequences: {total_seq_count}\n')
        f.write(f'# Promoters without custom criteria: {no_pattern_count}\n\n')
        f.write(f'# Share without custom criteria: {no_pattern_count/total_seq_count:.4f}\n\n')
        print("File successfully written:", output_filename)

def create_directory_if_not_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f'Directory "{dir_path}" created.')
    else:
        print(f'Directory "{dir_path}" already exists.')

if __name__ == "__main__":
    create_directory_if_not_exists(Filter_plot_directory_path)

    full_filepath = os.path.join(repeats_path, file_name)

    fn_list = file_name.split('_')
    output_file = f"{fn_list[0]}_{fn_list[1]}_prohib_nuc.txt"
    output_file = os.path.join(Filter_plot_directory_path, output_file)

    total_sequences, no_pattern_count = process_fasta_file()
    save_results_to_file(output_file, total_sequences, no_pattern_count)

