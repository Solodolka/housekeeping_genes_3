import clr
import os
import sys
import numpy as np
import re


# Указываем путь к папке с нашей DLL (текущая директория)
dll_path = os.path.abspath(".")
sys.path.append(dll_path)
# Загружаем библиотеку по имени файла (без расширения .dll)
clr.AddReference("DNAMotif_01")
# Импортируем класс из пространства имен C#
from DNAMotif_01 import DNAMotif, Motif

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


pfm_data = """
>MA0079.3	SP1
A [	2379	1323	2055	0	0	2751	73	0	1932	1877	1820	]
C [	734	0	0	0	0	4624	0	0	0	642	4271	]
G [	4969	6357	6679	8734	8734	0	8661	8734	6703	6215	1786	]
T [	652	1054	0	0	0	1359	0	0	99	0	857	]
"""
record = """
TCGCCATGAACTTGAGCCAGCCCAAGCCCTCCTGCAGCCAGCATGCTCTCCCACTCACTA
GTCTGACATCTGGCTACAAGGGGTGCTTGAGCCAAGAAAGATAAGTGCAGACTCCCGGAG
GGGACCAAGATTGCTGAAAGAGTTAAAAGGTTCTTTCTCCAAACCTTGGCCCTCCCTCTCCGGGGCGGGGC
CTGTGTCTACTTTATGGGTTTAATGAGAAACCCAGGCACCGCACAAGAGACTGCAGAAAT
AACTTCTGGGGCATCGCTAACTTATTCAGCAGCAAGAGACTTCAGGCTGCCAATGCTCCAGGGGGCGGGGC
GGGAAGATAACCCAGTGGAAAGAAGGTGGTTTTTACAGCCTGGGTCGGACGGAATTTTGC
AGAGCTGAGTGTGAGCTCCGCCCCTTACAGGACTGAAAAAGCCAGCTTAGGCTCCTCAGCGGGGGCGGGGC
CCCTCTTTGAGATATTTGATTAGGAATTTGGGGGCGGGCCCTGGCCTGGCACAGACTTGC
ACACTCTCCGTAGGCGCTCACTCCTCTCTTCCCTCTCTGCCACATTCTCCAACCCAAGGA
GACCAGACAG
"""

pwm_result = calculate_pwm_from_text(pfm_data)
print(pwm_result)
pwm_lists = pwm_result.flatten().tolist()

# Создаем экземпляр и вызываем методы
process_motif = DNAMotif()
process_motif.SetMatrix(pwm_lists)
process_motif.SetThreshold(0.0)
motifs = process_motif.FindMotifs(record)

for m in motifs:
    # Замените 'Sequence' и 'Position' на реальные имена полей вашей структуры Motif
    print(f"Найдено: {m.Weight} на позиции {m.Position}")
