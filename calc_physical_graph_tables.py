import os
from Bio import SeqIO

nuc_num = 2
tables_path = 'Physical_tables'
repeats_path = 'HK_promoters'
file_name = 'filtered_barley_epdnew_hv_promoters_all.fasta'
results_path = 'Phisycs'
table = {}  # Будет хранить таблицу
NUCLEOTIDES_IN_STRING = 90
aver_beg_values = [0] * NUCLEOTIDES_IN_STRING
aver_end_values = [0] * NUCLEOTIDES_IN_STRING
run_number = 1

NUCS = ['A', 'C', 'G', 'T']  # Список нуклеотидов для формирования ключа

def load_all_tables():
    global table, run_number
    # Получаем список всех файлов в папке с таблицами
    files = os.listdir(tables_path)

    for filename in files:
        table = {}
        # Формируем полный путь к файлу
        filepath = os.path.join(tables_path, filename)

        # Заполняем таблицу
        with open(filepath, 'r', encoding='utf-8') as f:
            n = 0
            for line in f:
                parts = line.strip().split('\t')
                for p in parts:
                    key = NUCS[n // 4] + NUCS[n % 4]
                    n += 1
                    value = float(p)

                    # Добавляем запись в соответствующую таблицу
                    table[key] = value
        prefix = filename.replace(".txt", "_")
        prefix = prefix + file_name.split('_')[1] + '_'
        run_number = 1

        read_and_process_repeats(repeats_path)
        write_results(aver_beg_values, prefix)
        print(aver_beg_values)


def read_and_process_repeats(directory):
    global run_number
    # Получаем список всех файлов в заданной директории
    filepath = os.path.join(directory, file_name)
    for record in SeqIO.parse(filepath, "fasta"):
        process_repeat_beg(record.seq)  # Обрабатываем начало
        run_number += 1


def process_repeat_beg(repeat):
    global NUCLEOTIDES_IN_STRING, run_number, aver_beg_values
    values = []
    key = repeat[:nuc_num]
    for i in range(0, shifts_num):
        values.append(table.get(key, 0))  # Берём значение из соответствующей таблицы
        key = key[1:] + repeat[nuc_num + i]
    # write_results(values, '_beg')

    # Вычисляем среднее значение
    coeff0 = (run_number - 1) / run_number
    coeff1 = 1 / run_number
    for i, v in enumerate(values, start=0):
        aver_beg_values[i] = coeff0 * aver_beg_values[i] + coeff1 * v


def write_results(values, pref=''):
    filepath = os.path.join(results_path, pref + str(nuc_num) + '.txt')
    with open(filepath, 'a', encoding='utf-8') as f:
        for v in values:
            v_str = f"{v:.8f}"  # Сохраняем число с точностью до 8 знаков после запятой
            f.write(v_str + '\t')
        f.write('\n')


load_all_tables()
