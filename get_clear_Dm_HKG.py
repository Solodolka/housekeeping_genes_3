import os

results_path = 'HK_genes'

def process_genes(input_file_path, output_file_path):
    # Чтение файла и получение списка уникальных генов
    with open(input_file_path, 'r') as input_file:
        genes = input_file.read().strip()

    unique_genes = sorted(set(genes.split(', ')))

    # Запись результата в новый файл
    with open(output_file_path, 'w') as output_file:
        for gene in unique_genes:
            output_file.write(f"{gene}\n")


input_file_path = os.path.join(results_path, 'Dm_HKG.txt')  # Исходный файл с генами
output_file_path = os.path.join(results_path, 'Dm_HKG_unique.txt')  # Выходной файл

process_genes(input_file_path, output_file_path)