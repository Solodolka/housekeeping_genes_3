import os
import matplotlib.pyplot as plt

# Директория с файлом
directory = 'Nucleotides_probabilities'
plot_directory = 'Nucleotides_probabilities_plots'
file_name = 'result_barley.txt'

def plot_nucleotides(file_path):
    """
    Функция строит графики вероятности нахождения четырех типов нуклеотидов,
    используя значения из указанного файла.

    Параметры:
        file_path (str): путь к файлу с данными
    """
    # Читаем данные из файла
    positions = []  # Список позиций (номер строки)
    probabilities = {'A': [], 'C': [], 'G': [], 'T': []}

    with open(file_path, 'r') as file:
        for i, line in enumerate(file):
            # Преобразуем строку в список чисел
            nums = list(map(float, line.strip().split('\t')))

            # Добавляем каждое число в соответствующий тип нуклеотида
            probabilities['A'].append(nums[0])
            probabilities['C'].append(nums[1])
            probabilities['T'].append(nums[2])
            probabilities['G'].append(nums[3])

            # Заполняем список позиций
            positions.append(i - 50)  # Индексация начинается с 1
            if i == 85:
                break

    # Строим графики
    plt.figure(figsize=(10, 6))
    colors = ['black', 'blue', 'yellow', 'red']
    labels = ['Аденин (A)', 'Цитозин (C)', 'Гуанин (G)', 'Тимин (T)']

    for nucleotide, color, label in zip(probabilities.keys(), colors, labels):
        plt.plot(positions, probabilities[nucleotide],
                 marker='.', markersize=0, color=color, label=label)

    # Установка собственных меток для оси X
    step = 10  # Метки ставятся через каждые 10 единиц
    custom_ticks = positions[::step]  # Берём каждую десятую позицию
    custom_labels = [str(i) if i != 0 else 'TSS' for i in custom_ticks]
    plt.xticks(custom_ticks, custom_labels)  # Применяем новые метки и позиции

    # Подписываем оси
    plt.axvline(x=0, color="red", linestyle="--", linewidth=0.5)  # Линия вертикальной разметки
    plt.axvline(x=-28, color="blue", linestyle="--", linewidth=0.5)  # Линия вертикальной разметки

    plt.legend(loc="upper right")
    plt.title('Вероятность наличия нуклеотидов по позициям')
    plt.xlabel('Позиция')
    plt.ylabel('Вероятность')
    plt.grid(True)

    # Сохраняем график
    output_filename = f'{os.path.splitext(os.path.basename(file_path))[0]}_graphs.png'
    output_filename = os.path.join(plot_directory, output_filename)
    plt.savefig(output_filename)
    print(f'Графики сохранены как {output_filename}')
    plt.close()


if __name__ == "__main__":
    # Полный путь к файлу
    full_path = os.path.join(directory, file_name)
    plot_nucleotides(full_path)
