# Fable, 15.09.2026. Загрузчик кода Bianchi–Padurariu (qc_g2_bielliptic.sage, GitHub bianchifrancesca/QC_bielliptic,
# скачан 15.09.2026 как текст) в обычный python3 + sage.all: препарсинг Sage-синтаксиса и exec в общее пространство имён.
import sys, os
from sage.all import *
from sage.repl.preparse import preparse_file
HERE=os.path.dirname(os.path.abspath(__file__))
src=open(os.path.join(HERE,'qc_g2_bielliptic_patched.sage')).read()
NS=globals()
exec(preparse_file(src), NS)
