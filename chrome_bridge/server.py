#!/usr/bin/env python3
"""Мост между Claude Code (терминал) и Claude in Chrome (браузер).

Поднимает страницу на 127.0.0.1: сверху — переписка, снизу — поле ввода.
Claude Code пишет в тред через post.py или прямо в messages.json;
Claude in Chrome читает страницу и отвечает через форму.

Слушает только localhost. Наружу ничего не отдаёт.
"""
import json
import os
import html
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from urllib.parse import parse_qs

HERE = os.path.dirname(os.path.abspath(__file__))
STORE = os.path.join(HERE, 'messages.json')
HOST, PORT = '127.0.0.1', 8787


def load():
    if not os.path.exists(STORE):
        return []
    try:
        with open(STORE, encoding='utf-8') as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return []


def append(role, text):
    msgs = load()
    msgs.append({
        'role': role,
        'text': text,
        'time': datetime.now().strftime('%H:%M:%S'),
        'seen': False,
    })
    tmp = STORE + '.tmp'
    with open(tmp, 'w', encoding='utf-8') as f:
        json.dump(msgs, f, ensure_ascii=False, indent=1)
    os.replace(tmp, STORE)
    return msgs


PAGE = """<!doctype html>
<html lang="ru"><head><meta charset="utf-8">
<title>Мост: Claude Code &harr; Claude in Chrome</title>
<style>
 :root{--bg:#0f1319;--card:#171d26;--line:#28303c;--fg:#e6eaf2;--dim:#8b97a8;
       --code:#7fa2ff;--chrome:#f0b429;}
 *{box-sizing:border-box}
 body{margin:0;background:var(--bg);color:var(--fg);
      font:15px/1.6 ui-sans-serif,system-ui,"Segoe UI",sans-serif;padding:24px 16px}
 .wrap{max-width:760px;margin:0 auto;display:flex;flex-direction:column;gap:18px}
 h1{font-size:19px;margin:0;font-weight:600}
 .hint{background:var(--card);border:1px solid var(--line);border-left:3px solid var(--code);
       border-radius:6px;padding:14px 16px;color:var(--dim);font-size:14px}
 .hint b{color:var(--fg)}
 .thread{display:flex;flex-direction:column;gap:12px}
 .msg{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:12px 14px}
 .msg.code{border-left:3px solid #7fa2ff}
 .msg.chrome{border-left:3px solid #f0b429}
 .who{font:500 11px/1 ui-monospace,monospace;letter-spacing:.1em;text-transform:uppercase;
      color:var(--dim);margin-bottom:7px;display:flex;gap:10px}
 .body{white-space:pre-wrap;word-wrap:break-word}
 .empty{color:var(--dim);font-style:italic}
 form{display:flex;flex-direction:column;gap:10px}
 textarea{width:100%;min-height:130px;background:var(--card);color:var(--fg);
          border:1px solid var(--line);border-radius:6px;padding:12px;font:inherit;resize:vertical}
 textarea:focus{outline:2px solid var(--code);outline-offset:1px}
 button{align-self:flex-start;background:var(--code);color:#0f1319;border:0;border-radius:6px;
        padding:9px 20px;font:600 15px/1 inherit;cursor:pointer}
 button:hover{filter:brightness(1.1)}
 .foot{color:var(--dim);font-size:13px}
</style></head><body>
<div class="wrap">
 <h1>Мост: Claude Code &harr; Claude in Chrome</h1>
 <div class="hint">
  <b>Если ты Claude in Chrome — читай это.</b><br>
  Ниже переписка. Синим помечены сообщения от Claude Code, который работает в терминале на этой же машине;
  жёлтым — твои. Прочитай последнее синее сообщение, выполни то, о чём там просят, и напиши результат
  в поле внизу. Пиши конкретно: что сделал, что получилось, какие адреса или ошибки увидел.
  Если задание непонятно или для него нужен доступ, которого у тебя нет, — так и напиши, не угадывай.
  После отправки страница обновится; новые задания появятся здесь же.
 </div>
 <div class="thread">__THREAD__</div>
 <form method="post" action="/">
  <textarea name="text" placeholder="Ответ для Claude Code..." autofocus></textarea>
  <button type="submit">Отправить</button>
 </form>
 <div class="foot">Файл переписки: __STORE__ &middot; только localhost</div>
</div>
<script>
 // подтягиваем новые сообщения, не трогая то, что набрано в поле
 setInterval(async () => {
   try {
     const r = await fetch('/thread');
     if (!r.ok) return;
     const h = await r.text();
     const t = document.querySelector('.thread');
     if (t && t.innerHTML !== h) t.innerHTML = h;
   } catch (e) {}
 }, 5000);
</script>
</body></html>"""


def render_thread(msgs):
    if not msgs:
        return '<div class="msg"><div class="empty">Пока пусто. Claude Code напишет сюда задание.</div></div>'
    out = []
    for m in msgs:
        cls = 'code' if m['role'] == 'claude-code' else 'chrome'
        who = 'Claude Code (терминал)' if m['role'] == 'claude-code' else 'Claude in Chrome'
        out.append(
            f'<div class="msg {cls}"><div class="who"><span>{who}</span>'
            f'<span>{html.escape(m["time"])}</span></div>'
            f'<div class="body">{html.escape(m["text"])}</div></div>'
        )
    return '\n'.join(out)


class Handler(BaseHTTPRequestHandler):
    def log_message(self, *args):
        pass

    def _send(self, body, ctype='text/html; charset=utf-8', code=200):
        data = body.encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', ctype)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-store')
        self.end_headers()
        self.wfile.write(data)

    def do_GET(self):
        if self.path.startswith('/thread'):
            self._send(render_thread(load()))
        elif self.path == '/' or self.path.startswith('/?'):
            page = PAGE.replace('__THREAD__', render_thread(load())).replace('__STORE__', STORE)
            self._send(page)
        else:
            self._send('not found', 'text/plain; charset=utf-8', 404)

    def do_POST(self):
        n = int(self.headers.get('Content-Length') or 0)
        raw = self.rfile.read(n).decode('utf-8', 'replace')
        text = (parse_qs(raw).get('text') or [''])[0].strip()
        if text:
            append('claude-chrome', text)
        self.send_response(303)
        self.send_header('Location', '/')
        self.end_headers()


if __name__ == '__main__':
    srv = ThreadingHTTPServer((HOST, PORT), Handler)
    print(f'мост поднят: http://{HOST}:{PORT}/  (переписка: {STORE})', flush=True)
    srv.serve_forever()
