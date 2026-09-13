# Fable, 14.09.2026. Параллельный прогон six_review.check_one по всем 1240 закрытиям этапа «шесть клеток».
import sys, json, time, multiprocessing as mp
sys.path.insert(0, '/home/kep/magicKube/fable_census_review')
import six_review as six

OUT = '/home/kep/magicKube/fable_census_review/six_review_par.jsonl'


def work(item):
    slope, abc, kind, meth = item
    t0 = time.time()
    try:
        rec = six.check_one(slope, abc, kind, tL=600, tM=120)
    except Exception as e:
        rec = {'slope': slope, 'abc': abc, 'kind': kind, 'closed': False, 'error': str(e)[:200]}
    rec['claude_method'] = meth
    rec['sec'] = round(time.time() - t0, 1)
    return rec


if __name__ == '__main__':
    rows = [json.loads(l) for l in open(six.SRC) if l.strip()]
    closed = [d for d in rows if d.get('closed')]
    done = set()
    for f in (six.OUT, OUT):
        try:
            for l in open(f):
                d = json.loads(l); done.add((d['slope'], tuple(d['abc']), d['kind']))
        except FileNotFoundError:
            pass
    todo = [(d['slope'], d['by'][0], d['by'][1], 'ellrank ' + json.dumps(d['info'])) for d in closed
            if (d['slope'], tuple(d['by'][0]), d['by'][1]) not in done]
    print('всего закрытий', len(closed), 'уже проверено', len(done), 'к проверке', len(todo), flush=True)
    n = 0
    with mp.get_context('fork').Pool(10) as pool:
        for rec in pool.imap_unordered(work, todo):
            with open(OUT, 'a') as f:
                f.write(json.dumps(rec, ensure_ascii=False) + '\n')
            n += 1
            if n % 20 == 0 or not rec['closed']:
                print(n, '/', len(todo), rec['slope'], rec['abc'], rec['kind'], 'closed=', rec['closed'],
                      'L=', rec.get('L_ratio'), 'mw=', rec.get('mwrank_bound'), flush=True)
    print('готово', flush=True)
