import subprocess
import os
import sys

# Список скриптов, которые нужно запустить
scripts_to_run = [
    # 'find_by_restrictions.py',
    'Calc_n_plot_all.py',
    # 'Calc_n_plot_all_test.py',
    # 'plot_cleavage_and_physics.py',
    # 'find_GC.py',
    # 'find_TATA.py',
]

animals = [
    'H_sapiens',
    'A_thaliana',
    'H_vulgare',
    'M_mulatta',
    'M_musculus',
    'Z_mays',
]

for animal in animals:
# Цикл по всем скриптам
    for script in scripts_to_run:
        # Создаем копию текущих переменных окружения и добавляем/обновляем нашу
        env = os.environ.copy()
        env['MY_GLOBAL_VAR'] = animal
        python_executable = sys.executable
        print(f"\nRun {script} with {animal}...")

        # Запускаем дочерний процесс с новым окружением
        # check=True вызовет ошибку, если скрипт завершится с ненулевым кодом
        subprocess.run([python_executable, script], env=env, check=True)