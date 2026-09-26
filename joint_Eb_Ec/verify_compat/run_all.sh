#!/bin/sh
# Воспроизведение независимой проверки контроля совместимости Codex (b=34, c=3400).
# Порядок: PARI ранги/корневое число/аналитический ранг; образующие Codex и насыщение;
# граница индекса (Sage, минимальная модель); формула Монски + Таннелл; доп. поиск пар.
set -e
cd "$(dirname "$0")"
gp -q ranks_pari.gp < /dev/null > ranks_pari.log 2>&1
gp -q sat_pari.gp   < /dev/null > sat_pari.log 2>&1
env DOT_SAGE=/tmp/claude_vc_sage python3 index_bound_sage.py > index_bound.log 2>&1
gp -q sel2_pari.gp  < /dev/null > sel2_pari.txt 2>&1
python3 monsky_tunnell.py > monsky_tunnell.log 2>&1
gp -q compat_search.gp < /dev/null > compat_search.log 2>&1
