def sort_lines(input_filename, output_filename):
    try:
        # Чтение строк из входного файла
        with open(input_filename, 'r', encoding='utf-8') as infile:
            lines = infile.readlines()

        # Удаление завершающих символов переноса строки (\n) и сортировка
        sorted_lines = sorted(line.strip() for line in lines)

        # Запись отсортированных строк в выходной файл
        with open(output_filename, 'w', encoding='utf-8') as outfile:
            outfile.write('\n'.join(sorted_lines))

        print("Файл успешно отсортирован.")
    except FileNotFoundError:
        print("Ошибка: указанный файл не найден.")
    except Exception as e:
        print(f"Произошла ошибка: {e}")


# Примеры использования
sort_lines('Clear_HK_genes.txt', 'Sorted_HK_genes.txt')