import os
from Bio import SeqIO
import matplotlib.pyplot as plt

repeats_path = os.getenv('MY_GLOBAL_VAR')
if repeats_path is None:
    repeats_path = 'Z_mays'

organizm = repeats_path.replace('_', '')

file_name_list = [
    f'all_{organizm}_epdnew_promoters.fasta',
    f'TATA_{organizm}_epdnew_promoters.fasta',
    f'TATAless_{organizm}_epdnew_promoters.fasta',
    # f'TATA_{organizm}_epdnew_promoters.fasta',
    # f'TATAless_{organizm}_epdnew_promoters.fasta',
]
NUCLEOTIDES_IN_STRING = 80  # Number of positions (nucleotides) to process

cleavage_tables_path = 'Cleavage_tables'
cleavage_txt_directory = 'Cleavage_txt'
Cleavage_txt_directory_path = repeats_path + "\\" + cleavage_txt_directory
cleavage_plot_directory = 'Cleavage_plots'
Cleavage_plot_directory_path = repeats_path + "\\" + cleavage_plot_directory

physical_tables_path = 'Physical_tables'
physical_txt_directory = 'Physical_txt'
Physicale_txt_directory_path = repeats_path + "\\" + physical_txt_directory
physical_plot_directory = 'Physical_plots'
Physicale_plot_directory_path = repeats_path + "\\" + physical_plot_directory

table = {}
aver_beg_values = []
run_number = 1
nuc_num = 2
NUCS = ['A', 'C', 'G', 'T']  # Список нуклеотидов для формирования ключа

def read_table():  # Загрузка таблицы 'AA...'
    filepath = os.path.join(cleavage_tables_path, 'nuc_table_' + str(nuc_num) + '.txt')
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            str_lst = line.split('\t')
            table[str_lst[0]] = float(str_lst[1])

def read_and_process_repeats():
    global run_number
    filepath = os.path.join(repeats_path, file_name)
    for record in SeqIO.parse(filepath, "fasta"):
        process_repeat(record.seq)  # Обрабатываем начало
        run_number += 1


def process_repeat(repeat):  # Обрабатываем начало
    global NUCLEOTIDES_IN_STRING, run_number, aver_beg_values
    values = []
    key = repeat[:nuc_num]
    for i in range(0, NUCLEOTIDES_IN_STRING - nuc_num + 1):
        values.append(table.get(key, 0))  # Если ключа нет, присваиваем 0
        key = key[1:] + repeat[nuc_num + i]

    # Расчёт среднего значения для каждой позиции
    coeff0 = (run_number - 1) / run_number
    coeff1 = 1 / run_number
    for i, v in enumerate(values, start=0):
        aver_beg_values[i] = coeff0 * aver_beg_values[i] + coeff1 * v

def txt_nucleotides(data, output_filename):
    with open(output_filename, 'a', encoding='utf-8') as f:
        for v in data:
            v_str = f"{v:.8f}"  # Сохраняем число с точностью до 8 знаков после запятой
            f.write(v_str + '\t')
        f.write('\n')
    print(f'Txt saved as {output_filename}')

def plot_nucleotides(data, output_filename, prefix):
    # Создание массива осей X с новыми индексами
    x_positions = list(range(-50, 25))  # [-50, -49, ..., 0 (TSS), ..., +35]

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

def process_all_physical_tables():
    global table, run_number, aver_beg_values, nuc_num

    nuc_num = 2
    # Получаем список всех файлов в папке с таблицами
    table_files = os.listdir(physical_tables_path)

    for table_filename in table_files:
        table = {}
        # Формируем полный путь к файлу
        filepath = os.path.join(physical_tables_path, table_filename)
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

#        print(table_filename)
#        print(table)

        run_number = 1
        aver_beg_values = [0] * (NUCLEOTIDES_IN_STRING - nuc_num + 1)
        read_and_process_repeats()

        prefix = table_filename.replace(".txt", "")
        output_filename = f"{prefix} " + fn_list[0] + '_' + fn_list[1] + ".png"
        output_filename = os.path.join(Physicale_plot_directory_path, output_filename)
        plot_nucleotides(aver_beg_values, output_filename, ". " + prefix)

        output_filename = f"{prefix} " + fn_list[0] + '_' + fn_list[1] + ".txt"
        output_filename = os.path.join(Physicale_txt_directory_path, output_filename)
        txt_nucleotides(aver_beg_values, output_filename)

def process_all_cleavage_tables():
    global run_number, aver_beg_values
    for i in range(3):
        nuc_num = (i + 1) * 2
        run_number = 1  # Возвращаем счётчик обратно к 1
        aver_beg_values = [0] * (NUCLEOTIDES_IN_STRING - nuc_num + 1)
        read_table()
        read_and_process_repeats()
        if (nuc_num == 6):
            prefix = '. DNase I cleavage'
        if (nuc_num == 4):
            prefix = '. Utrasonic 4 cleavage'
        if (nuc_num == 2):
            prefix = '. Utrasonic 2 cleavage'
        plot_nucleotides(aver_beg_values, prefix)

def create_directory_if_not_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f'Directory "{dir_path}" created.')
    else:
        print(f'Directory "{dir_path}" already exists.')


if __name__ == "__main__":
    for file_name in file_name_list:
        fn_list = file_name.split('_')
        create_directory_if_not_exists(Cleavage_txt_directory_path)
        create_directory_if_not_exists(Cleavage_plot_directory_path)
        for i in range(3):
            nuc_num = (i + 1) * 2
            run_number = 1  # Возвращаем счётчик обратно к 1
            aver_beg_values = [0] * (NUCLEOTIDES_IN_STRING - nuc_num + 1)
            read_table()
            read_and_process_repeats()
            if (nuc_num == 6):
                prefix = '. DNase I cleavage'
            if (nuc_num == 4):
                prefix = '. Utrasonic 4 cleavage'
            if (nuc_num == 2):
                prefix = '. Utrasonic 2 cleavage'

            output_filename = f"Cleavage {nuc_num} " + fn_list[0] + '_' + fn_list[1] + ".png"
            output_filename = os.path.join(Cleavage_plot_directory_path, output_filename)
            plot_nucleotides(aver_beg_values, output_filename, prefix)

            output_filename = f"Cleavage {nuc_num} " + fn_list[0] + '_' + fn_list[1] + ".txt"
            output_filename = os.path.join(Cleavage_txt_directory_path, output_filename)
            txt_nucleotides(aver_beg_values, output_filename)

        create_directory_if_not_exists(Physicale_txt_directory_path)
        create_directory_if_not_exists(Physicale_plot_directory_path)
        process_all_physical_tables()