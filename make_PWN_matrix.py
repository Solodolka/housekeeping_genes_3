import numpy as np
import re


def calculate_pwm_from_text(pfm_text, background=0.25, pseudocount=0.1):
    # Парсит текстовую PFM и преобразует её в PWM (Log-Odds).
    # 1. Извлекаем строки данных (A, C, G, T)
    # Ищем символ нуклеотида и всё, что внутри скобок [ ... ]
    rows = re.findall(r'([ACGT])\s+\[\s*([\d\s]+)\]', pfm_text)
    # Создаем словарь, чтобы гарантировать порядок A, C, G, T
    pfm_dict = {char: [float(x) for x in values.split()] for char, values in rows}
    # Превращаем в 2D массив в строгом порядке ACGT
    pfm = np.array([pfm_dict[n] for n in "ACGT"])
    # 2. Добавляем псевдоотсчеты
    pfm_adj = pfm + pseudocount
    # 3. Превращаем в вероятности (PPM)
    column_sums = pfm_adj.sum(axis=0)
    ppm = pfm_adj / column_sums
    # 4. Вычисляем Log-Odds (PWM)
    pwm = np.log2(ppm / background)
    return pwm


# Данные из прошлого запроса
pfm_data = """
>MA0079.3	SP1
A  [   857      0     99      0      0   1359      0      0      0   1054    652 ]
C  [  1786   6215   6703   8734   8661      0   8734   8734   6679   6357   4969 ]
G  [  4271    642      0      0      0   4624      0      0      0      0    734 ]
T  [  1820   1877   1932      0     73   2751      0      0   2055   1323   2379 ]
"""

pwm_result = calculate_pwm_from_text(pfm_data)

print("Матрица PWM (Log-Odds) для SP1:")
print(np.round(pwm_result, 2))
