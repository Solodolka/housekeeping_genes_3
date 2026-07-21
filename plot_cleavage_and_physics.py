import os
import matplotlib.pyplot as plt

repeats_path = os.getenv('MY_GLOBAL_VAR')
if repeats_path is None:
    repeats_path = 'Z_mays'

NUCLEOTIDES_IN_STRING = 80  # Number of positions (nucleotides) to process

cleavage_txt_directory = 'Cleavage_txt'
Cleavage_txt_directory_path = repeats_path + "\\" + cleavage_txt_directory
cleavage_plot_directory = 'Comparision_plots'
Cleavage_plot_directory_path = repeats_path + "\\" + cleavage_plot_directory

physical_txt_directory = 'Physical_txt'
Physicale_txt_directory_path = repeats_path + "\\" + physical_txt_directory
physical_plot_directory = 'Comparision_plots'
Physicale_plot_directory_path = repeats_path + "\\" + physical_plot_directory


def read_file_data(file_path):
    """Читает данные из файла и возвращает список чисел."""
    with open(file_path, 'r') as file:
        data = []
        for line in file:
            numbers = list(map(float, line.strip().split('\t')))
            data.extend(numbers)
    return data


def format_filename(filename):
    """Форматирует имя файла для подписи на графике."""
    name = filename.replace('.txt', '')
    name = name.replace('_2_', ' 2 ')
    name = name.replace('_4_', ' 4 ')
    name = name.replace('_6_', ' 6 ')
    name = name.replace('_2', ' ')
    name = name.replace('_', ' ')
    return name


def plot_files(file_list, table_name, path):
    x_positions = list(range(-50, 25))  # [-50, -49, ..., 0 (TSS), ..., +25]

    plt.figure(figsize=(12, 6))

    for file_postfix in file_list:
        file_name = table_name + ' ' + file_postfix + '.txt'
        file_name = os.path.join(path, file_name)
        data = read_file_data(file_name)
        plot_data = data[:len(x_positions)]
        plt.plot(x_positions, plot_data, marker='.', markersize=4, label=file_postfix)

    # Настройка осей и внешнего вида
    step = 10
    custom_ticks = x_positions[::step]
    custom_labels = [str(i) if i != 0 else 'TSS' for i in custom_ticks]
    plt.xticks(custom_ticks, custom_labels)

    plt.axvline(x=0, color="red", linestyle="--", linewidth=0.8, label='TSS')
    plt.axvline(x=-28, color="blue", linestyle="--", linewidth=0.8, label='-28')

    plt.xlabel('Relative Position to TSS')
    plt.ylabel('Value')
    if table_name == 'Cleavage 6':
        table_name = 'DNase I cleavage'
    plt.title(table_name + ' comparison')
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

    output_path = table_name + ' Comparison.png'
    output_path = os.path.join(Physicale_plot_directory_path, output_path)
    plt.savefig(output_path, bbox_inches='tight')
    print(f'Graph is saved as {output_path}')
    plt.close()


filter_list = [
    'all',
    'TATA',
    'TATAless',
    # 'TATA',
    # 'TATAless',
    # Добавьте другие файлы по необходимости
]

def create_directory_if_not_exists(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        print(f'Directory "{dir_path}" created.')
    else:
        print(f'Directory "{dir_path}" already exists.')


create_directory_if_not_exists(Cleavage_plot_directory_path)
create_directory_if_not_exists(Physicale_plot_directory_path)

physical_tables_path = 'Physical_tables'
table_files = os.listdir(physical_tables_path)
for i in range (len(filter_list)):
    filter_list[i] = filter_list[i] + '_' + repeats_path.replace('_', '')

for table_filename in table_files:
    table_name = table_filename.replace(".txt", "")
    plot_files(filter_list, table_name, Physicale_txt_directory_path)

for i in range(3):
    nuc_num = (i + 1) * 2
    table_name = f"Cleavage {nuc_num}"
    plot_files(filter_list, table_name, Cleavage_txt_directory_path)
