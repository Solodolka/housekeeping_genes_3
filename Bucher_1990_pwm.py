import numpy as np

# Данные частот (из вашей таблицы)
counts = {
    'A': [61, 16, 352, 3, 354, 268, 360, 222, 155, 56, 83, 82, 82, 68, 77],
    'C': [145, 46, 0, 10, 0, 0, 3, 2, 44, 135, 147, 127, 118, 107, 101],
    'G': [152, 18, 2, 2, 5, 0, 20, 44, 157, 150, 128, 128, 128, 139, 140],
    'T': [31, 309, 35, 374, 30, 121, 6, 121, 33, 48, 31, 52, 61, 75, 71]
}

def calculate_pwm(counts, pseudo=1):
    bases = 'ACGT'
    data = np.array([counts[b] for b in bases])

    # Сумма по столбцам (общее кол-во последовательностей)
    col_sums = data.sum(axis=0)

    # Вероятности с псевдоотсчетами
    # p = (n + s/4) / (N + s), где s - псевдоотсчет
    probs = (data + pseudo / 4) / (col_sums + pseudo)

    # Фоновая частота (0.25 для генома)
    background = 0.25

    # W = log2(P_observed / P_background)
    # В таблице Бухера значения часто нормализованы (max = 0)
    weights = np.log2(probs / background)

    # Нормализация: вычитаем максимум в каждом столбце (как в примере)
    weights -= weights.max(axis=0)

    return {base: weights[i].round(2) for i, base in enumerate(bases)}


def calculate_score(sequence, pwm):
    # Последовательность должна быть той же длины, что и PWM (15 символов)
    score = 0
    for i, base in enumerate(sequence.upper()):
        if base in pwm:
            score += pwm[base][i]
        else:
            score += -10 # Штраф за неопределенный нуклеотид (N)
    return score

# Пример (TATA-box):
seq = "GCTATAAAAGGGGGG"
pwm = calculate_pwm(counts)
positions = list(range(-3, 12))
print(f"{'Base':<5} | " + " | ".join(f"{p:^6}" for p in positions))
print("-" * (7 + 9 * len(positions)))
# Данные
for base in "ACGT":
    row_values = " | ".join(f"{v: >6.2f}" for v in pwm[base])
    print(f"{base:<5} | {row_values}")


current_score = calculate_score(seq, pwm)
print(f"Score: {current_score:.2f}")
print("Promoter detected:" if current_score >= -8.16 else "Not a promoter")
