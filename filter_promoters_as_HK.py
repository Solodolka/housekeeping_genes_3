from Bio import SeqIO


# Чтение списка генов из первого файла
def read_gene_list(file_path):
    genes = set()
    with open(file_path, 'r') as f:
        for line in f:
            # Ген расположен перед первым пробелом, приводим к верхнему регистру
            gene_name = line.split()[0].upper().strip()
            genes.add(gene_name)
    return genes


# Функция фильтрации последовательностей
def filter_sequences(fasta_file, output_file, genes):
    records_to_write = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Обрабатываем описание последовательности
        desc_parts = record.description.strip().split(' ')
        if len(desc_parts) >= 2:
            # Берём второе поле (название гена) и сравниваем с нашими генами
            gene_in_fasta = desc_parts[1].upper().strip()
            # Проверяем, начинается ли ген из Fasta с любого из наших генов
            for gene in genes:
                if gene_in_fasta.startswith(gene):
                    records_to_write.append(record)
                    break  # останавливаемся, если нашли хотя бы одно совпадение

    # Записываем отобранные последовательности в выходную фаил-фасту
    SeqIO.write(records_to_write, output_file, "fasta")


if __name__ == "__main__":
    input_gene_dir = "HK_genes"
    input_gene_file = "Dm_HKG_unique.txt"
    fasta_input_dir = "EPDnew_promoters"
    fasta_input_file = "drosophila_epdnew_dm6_promoters_all.fasta"
    output_dir = "HK_promoters"
    output_file = "Dm_HKG_" + fasta_input_file

    # Читаем список генов
    genes = read_gene_list(input_gene_dir + "\\" + input_gene_file)

    # Фильтруем записи и записываем их в выходной файл
    filter_sequences(fasta_input_dir + "\\" + fasta_input_file,
                     output_dir + "\\" + output_file,
                     genes)