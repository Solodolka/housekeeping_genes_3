from Bio import SeqIO

# Исходный файл и выходной файл
directory = 'M_musculus\\'
input_file = 'CommonHK_Mmusculus_epdnew_promoters.fasta'
output_file = 'trimmed_CommonHK_Mmusculus_epdnew_promoters.fasta'


def contains_invalid_chars(seq):
    """Проверяет, содержит ли последовательность недопустимые символы"""
    valid_chars = {'A', 'T', 'G', 'C', 'N'}
    return any(char.upper() not in valid_chars for char in seq)


def trim_sequence(input_fasta, output_fasta):
    # Чтение файла FASTA
    records = list(SeqIO.parse(directory + input_fasta, "fasta"))
    print(f"Длина записи:  {len(records[0])}")


    # Открываем файл для записи сразу
    with open(output_fasta, "w") as out_handle:
        for record in records:
            if len(record.seq) != len(records[0]):
                print(f"Пропущена запись {record.id}: длина {len(record.seq)}.")
                continue

            # Проверяем, содержит ли последовательность недопустимые символы
            if contains_invalid_chars(record.seq):
                print(f"Пропущена запись {record.id}: содержит недопустимые символы.")
                continue

            # Обрезаем нужный участок с -50 по +30
            trimmed_seq = record.seq[0:81]
            # Формируем новое описание
            new_description = f"{record.description.split(';')[0]} ; range -50 to 30"
            # Создаем новый объект записи
            new_record = record[:]
            new_record.seq = trimmed_seq
            new_record.description = new_description

            # Записываем одну запись в открытый файл
            SeqIO.write([new_record], out_handle, "fasta")


if __name__ == "__main__":
    trim_sequence(input_file, directory + output_file)
