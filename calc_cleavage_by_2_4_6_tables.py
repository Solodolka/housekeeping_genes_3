import os
from Bio import SeqIO

nuc_num = 6
tables_path = 'Cleavage_tables'
repeats_path = 'HK_promoters'
file_name = 'filtered_barley_epdnew_hv_promoters_all.fasta'
results_path = 'Cleavage'
table = {}
NUCLEOTIDES_IN_STRING = 95
aver_beg_values = [0] * NUCLEOTIDES_IN_STRING
run_number = 1


def read_table():  #Загрузка таблицы 'AA...' - 3.14
    filepath = os.path.join(tables_path, 'nuc_table_' + str(nuc_num) + '.txt')
    with open(filepath, 'r', encoding='utf-8') as f:
        # Читаем файл построчно
        for line in f:
            str_lst = line.split('\t')
            table[str_lst[0]] = float(str_lst[1])

def read_and_process_repeats(directory): # Перебор всех файлов с репитами
    global run_number
    # Получаем список всех файлов в заданной директории
    filepath = os.path.join(directory, file_name)
    for record in SeqIO.parse(filepath, "fasta"):
        process_repeat_beg(record.seq)  # Обрабатываем начало
        run_number += 1

def process_repeat_beg(repeat):   # Обрабатываем начало
    global NUCLEOTIDES_IN_STRING, run_number, aver_beg_values
    values = []
    key = repeat[:nuc_num]
    for i in range(0, shifts_num):
        values.append(table[key])
        key = key[1:] + repeat[nuc_num + i]
    # write_results(values)

    # calc average value for each position
    coeff0 = (run_number-1)/run_number
    coeff1 = 1/run_number
    for i, v in enumerate(values, start=0):
        aver_beg_values[i] = coeff0 * aver_beg_values[i] + coeff1 * v
    print (values)

def write_results(values):
    prefix = '_' + file_name.split('_')[1]
    filepath = os.path.join(results_path, 'Ultrasonic cleavage_' + str(nuc_num) + prefix + '.txt')
    with open(filepath, 'a', encoding='utf-8') as f:
        for v in values:
            v_str = f"{v:.8f}"
            # v_str = v_str.replace('.', ',') # FOR RUSSIAN Excel
            f.write(v_str + '\t')
        f.write('\n')


read_table()
read_and_process_repeats(repeats_path)
write_results(aver_beg_values)
print(aver_beg_values)
