
import sys, json
sys.path.insert(0, "/home/kep/magicKube/isogeny_mult/fable")
from rank_classes import work
for sl in sys.argv[2:]:
    r, s = map(int, sl.split("/")); res = work((r, s))
    with open(sys.argv[1], "a") as f: f.write(json.dumps(res) + "\n")
