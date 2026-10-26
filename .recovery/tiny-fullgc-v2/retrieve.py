import hashlib
import json
import shutil
import sys
from pathlib import Path

scratch, destination = map(Path, sys.argv[1:3])
destination.mkdir(parents=True, exist_ok=False)
selected = list(scratch.glob('*.log')) + list(scratch.glob('*.txt'))
selected += [scratch / 'inputs.sha256', scratch / 'CMakeLists.txt', scratch / 'src']
if (scratch / 'build/Testing').is_dir():
    selected.append(scratch / 'build/Testing')
records = {}
for source in selected:
    if not source.exists():
        continue
    relative = source.relative_to(scratch)
    target = destination / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(source, target)
        paths = [p for p in source.rglob('*') if p.is_file()]
    else:
        shutil.copy2(source, target)
        paths = [source]
    for path in paths:
        relative = path.relative_to(scratch)
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != hashlib.sha256((destination / relative).read_bytes()).hexdigest():
            raise RuntimeError('Retrieval checksum mismatch: ' + str(relative))
        records[str(relative)] = digest
(destination / 'retrieval-checksums.json').write_text(json.dumps(records, indent=2) + '\n')
print('VERIFIED_RETRIEVAL ' + str(destination), flush=True)
shutil.rmtree(scratch)
