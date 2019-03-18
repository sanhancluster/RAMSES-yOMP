import sh
from glob import glob
from tqdm import tqdm

files = glob('**/*.f90')

FROM = [r'\.gt\.', r'\.lt\.', r'\.ne\.', r'\.le\.', r'\.ge\.']
TO   = ['>',      '<',      '\/=',    '<=',     '>=']
for _from, _to in zip(FROM, TO):
    for f in tqdm(files):
        arg = (r'-i -e"s/%(fr)s/%(to)s/g"' %
               dict(fr=_from, to=_to))
        print(arg.split())
        sh.sed(arg.split() + [f])
