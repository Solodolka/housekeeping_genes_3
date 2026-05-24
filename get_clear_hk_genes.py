def extract_and_save_info(input_file_path, output_file_path):
    with open(input_file_path, 'r', encoding='utf-8') as input_file, \
            open(output_file_path, 'w', encoding='utf-8') as output_file:

        for line in input_file:
            # Разделение строки по двум элементам '**'
            parts = line.strip().split('**')

            if len(parts) >= 3:
                # Вторая часть списка — это нужная нам информация
                extracted_data = parts[1]

                # Записываем извлечённую информацию в выходной файл
                output_file.write(f"{extracted_data}\n")


# Пример использования
input_file_name = 'HK_proteins_description.txt'  # исходный файл
output_file_name = 'Clear_HK_genes.txt'  # результирующий файл
extract_and_save_info(input_file_name, output_file_name)