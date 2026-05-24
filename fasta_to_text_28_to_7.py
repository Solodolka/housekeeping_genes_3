from Bio import SeqIO

# Открываем исходный FASTA-файл для чтения
input_file = 'HK_promoters//Dm_HKGdrosophila_epdnew_dm6_promoters_all.fasta'
output_file = 'trimmed_sequences.txt'

with open(input_file, 'r') as handle, open(output_file, 'w') as out_handle:
    # Чтение записей из FASTA-файла
    for record in SeqIO.parse(handle, 'fasta'):
        # Разбиваем описание на элементы и получаем имя гена
        gene_name = record.description.split()[1]

        # Вырезаем участок последовательности от -28 до -7
        trimmed_seq = str(record.seq)[-28:-7]

        # Запись имени гена и обрезанной последовательности в выходной файл
        out_handle.write(f"{gene_name}\n")
        out_handle.write(trimmed_seq + "\n\n")
