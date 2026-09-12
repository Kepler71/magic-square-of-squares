#!/usr/bin/env python3
"""Отправить сообщение браузерному агенту: python3 post.py "текст" (или через stdin)."""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from server import append, STORE  # noqa: E402

text = ' '.join(sys.argv[1:]).strip() or sys.stdin.read().strip()
if not text:
    sys.exit('пусто — нечего отправлять')
append('claude-code', text)
print(f'отправлено ({len(text)} символов) → {STORE}')
