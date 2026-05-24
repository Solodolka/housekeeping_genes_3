import os
from collections import Counter
from Bio import SeqIO

# Переменные и пути
repeats_path = 'HK_promoters'  # Каталог с файлами FASTA
file_name = 'filtered_drosophila_epdnew_dm6_promoters_all.fasta'  # Имя файла FASTA
results_path = 'Nucleotides_probabilities'  # Каталог для вывода результатов
output_file = os.path.join(results_path, 'kmer_frequencies.txt')  # Название выходного файла


def extract_kmers(seq):
    """
    Извлекаем последовательность из заданного диапазона (позиции 20-27).
    Возвращает строку длины 8 нуклеотидов.
    """
    start_pos = 20
    end_pos = 27
    kmer = seq[start_pos:end_pos + 1].upper()  # Используем верхний регистр для единообразия
    return str(kmer)


def main():
    kmers_counter = Counter()

    # Читаем файл и извлекаем последовательности длиной 8 нуклеотидов
    fasta_file = os.path.join(repeats_path, file_name)
    for record in SeqIO.parse(fasta_file, "fasta"):
        kmer = extract_kmers(record.seq)
        if len(kmer) != 8: continue  # Проверяем длину строки
        kmers_counter[kmer] += 1  # Увеличиваем счётчик найденной комбинации

    # Сортируем полученные комбинации по убыванию частоты
    sorted_kmers = sorted(kmers_counter.items(), key=lambda x: x[1], reverse=True)

    # Записываем результат в файл
    with open(output_file, 'w', encoding='utf-8') as output:
        for kmer, freq in sorted_kmers:
            output.write(f'{kmer}\t{freq}\n')


if __name__ == "__main__":
    main()