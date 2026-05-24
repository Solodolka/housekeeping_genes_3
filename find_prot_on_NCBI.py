from Bio import Entrez

# Устанавливаем электронную почту для связи с NCBI
Entrez.email = "soldatenkova.o@gmail.com"

# Список символов белков
symbols = ["Smim10l2a_1", "Adam17_1", "DYNC2LI1_1", "TUBB3_1"]

# Создаем пустые списки для хранения результатов
names = []
id = []

# Ищем каждую запись по символу белка
for symbol in symbols:
    # Осуществляем поиск по базе данных Proteins
    handle = Entrez.esearch(db="protein", term=symbol)
    record = Entrez.read(handle)
    id_list = record['IdList']

    if len(id_list) > 0:
        # Получаем детальную информацию по каждому найденному идентификатору
        handle = Entrez.esummary(db="protein", id=id_list[0])  # Выбираем первый результат
        summary_record = Entrez.read(handle)

        names.append(summary_record[0]['Title'])
        id.append(summary_record[0]['Id'])
    else:
        names.append("Не найден")
        id.append("")

# Выводим результаты
for i, symbol in enumerate(symbols):
    print(f"{symbol}_1:  Id - {id[i]}, Название - {names[i]}")