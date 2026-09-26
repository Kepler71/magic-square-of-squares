from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
tasks=[('lift','class_lists.py'),('lift','check_lift_support.py'),
       ('index','index_checks.py'),('index','certify_local_relation_positive.py'),
       ('claude_review','verify_independence_certificate.py'),
       ('claude_review','bounded_local_review.py')]
for folder,script in tasks:
    print('Checking',folder,script,flush=True)
    subprocess.run([sys.executable,str(root/folder/script)],cwd=root/folder,check=True)
print('All listed controls passed. This does not assert global nonexistence.')
