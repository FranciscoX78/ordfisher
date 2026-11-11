import subprocess,sys,os
from pathlib import Path
import json

def test_logistic_small(tmp_path):
    root=Path(__file__).resolve().parents[1]
    env=dict(os.environ);env['PYTHONPATH']=str(root)+os.pathsep+env.get('PYTHONPATH','')
    r=subprocess.run([sys.executable,str(root/'scripts/run_logistic.py'),
              '--out',str(tmp_path/'logistic'),'--seed','123','--repetitions','2',
              '--retained','300','--burn','200',
              '--r-min','3.70','--r-max','3.71','--r-step','0.01'],
              capture_output=True,text=True,cwd=root,env=env)
    assert r.returncode==0,r.stderr
    assert (tmp_path/'logistic/logistic_summary.csv').exists()
    assert (tmp_path/'logistic/logistic_metadata.json').exists()


def test_verify_markov_script(tmp_path):
    root=Path(__file__).resolve().parents[1]
    path=tmp_path/'report.json'
    env=dict(os.environ);env['PYTHONPATH']=str(root)+os.pathsep+env.get('PYTHONPATH','')
    r=subprocess.run([sys.executable,str(root/'scripts/verify_markov.py'),'--out',str(path)],
                     capture_output=True,text=True,cwd=root,env=env)
    assert r.returncode==0,r.stderr
    j=json.loads(path.read_text())
    assert abs(j['spectral_numerical']['tv_ofs_normalized']-.23122928)<1e-7
