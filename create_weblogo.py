from weblogo import *

# 1. Загрузка последовательностей из FASTA файла
file_directory = "HK_promoters"
file_name = "trimmed_barley_epdnew_hv_promoters_all.fasta"
output_filename = os.path.join(file_directory, file_name)
fin = open(output_filename, "r")
seqs = read_seq_data(fin)


# 2. Создание логотипа
logo = LogoData.from_seqs(seqs)
options = LogoOptions()
options.logo_title = "Barley WebLogo"
# Настройка параметров логотипа
# делаем столько стопок, сколько нужно, чтобы занять одну строку
options.stacks_per_line = 81
options.stack_width = std_sizes['medium']
options.stack_aspect_ratio = 30
options.color_scheme = ColorScheme([
    SymbolColor("G", "yellow"),  # Guanine (Гуанин) зеленым
    SymbolColor("C", "blue"), # Cytosine (Цитозин) желтым
    SymbolColor("A", "black"),    # Adenine (Аденин) красным
    SymbolColor("T", "red")    # Thymine (Тимин) синим
])
options.yaxis_scale = 0.4
options.fineprint_font_size = 30
options.logo_height = 5
options.resolution = 600
options.unit_name = "bits"
options.stack_spacing = 0              # Минимизирует расстояние между стопками
options.show_fineprint = True         # Отключает подпись снизу
options.logo_margin_top = 0            # Убирает верхний отступ
options.logo_margin_bottom = 0         # Убирает нижний отступ
format = LogoFormat(logo, options)

# 3. Сохранение в файл (формат PDF по умолчанию)
image_bytes = jpeg_formatter(logo, format)
with open("my_logo.jpg", "wb") as f:
    f.write(image_bytes)
