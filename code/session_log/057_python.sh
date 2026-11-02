# session cell 57 | 2026-10-15 21:29 | env: python | language: bash
# files written: []

uname -m; sw_vers 2>/dev/null | head -3; df -h . | tail -1
for v in 2faf35be-0a80-4e7b-9433-1cc8a4077080 61fdc6b1-6d5b-4a7b-844a-1d2d354aba93 1d0f969c-960e-48fd-b8ab-e35e28955459 ab02abc4-bcb5-4de5-8b24-d88541cb0e21 75805d05-540f-4d38-b3b2-258490811841 cf52646c-30dc-48dd-baa0-2812086eee42 24f22222-c59f-4660-a6c9-b1943330d6f0; do :; done
python - <<'E'
import glob
E
curl -s -m 20 "https://conda.anaconda.org/bioconda/osx-arm64/repodata.json" -o /tmp/arm.json; python -c "
import json;d=json.load(open('/tmp/arm.json'));p={**d.get('packages',{}),**d.get('packages.conda',{})};print('arm64 antismash',sorted({v['version'] for v in p.values() if v['name']=='antismash'})[-3:])"
curl -s -m 60 "https://conda.anaconda.org/bioconda/osx-64/repodata.json" -o /tmp/x64.json; python -c "
import json;d=json.load(open('/tmp/x64.json'));p={**d.get('packages',{}),**d.get('packages.conda',{})};print('osx-64 antismash',sorted({v['version'] for v in p.values() if v['name']=='antismash'})[-3:])"
/usr/bin/arch -x86_64 /usr/bin/true && echo rosetta_ok