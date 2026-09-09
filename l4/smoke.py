# Experimental AENET runtime check; MPL-2.0, see src/license-header.txt.
import math, os, pathlib, re, subprocess, sys
root=pathlib.Path(sys.argv[1]).resolve()
work=pathlib.Path(sys.argv[2]).resolve();work.mkdir(parents=True,exist_ok=True)
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
def run(name,*args):
 exe=next((root/'bin').glob(name+'.x*'))
 r=subprocess.run([str(exe),*args],cwd=work,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=90)
 (work/(name+'.log')).write_text(r.stdout)
 assert r.returncode==0 and 'Error:' not in r.stdout,(name,r.stdout[-2000:])
 return r.stdout
(work/'Cu.stp').write_text('DESCR\nSmoke\nEND DESCR\nATOM Cu\nENV 1\nCu\nRMIN 0.5\nBASIS type=chebyshev\nradial_Rc=5 radial_N=3 angular_Rc=5 angular_N=2\n')
for i in range(8):
 d=2.+i*.15;e=(d-2.5)**2-1
 (work/f'{i}.xsf').write_text(f'# total energy = {e} eV\nATOMS\nCu 0 0 0 {2*(d-2.5)} 0 0\nCu {d} 0 0 {-2*(d-2.5)} 0 0\n')
(work/'generate.in').write_text('OUTPUT smoke.train\nTYPES\n1\nCu 0\nSETUPS\nCu Cu.stp\nFILES\n8\n'+''.join(f'{i}.xsf\n' for i in range(8)))
run('generate','generate.in');assert (work/'smoke.train').stat().st_size>0
(work/'train.in').write_text('TRAININGSET smoke.train\nTESTPERCENT 0\nITERATIONS 2\nMETHOD\nbfgs\nNETWORKS\nCu Cu.nn 1 3:tanh\n')
run('train','train.in');assert (work/'Cu.nn').stat().st_size>0
(work/'predict.in').write_text('TYPES\n1\nCu\nNETWORKS\nCu Cu.nn\nFORCES\n')
out=run('predict','predict.in','3.xsf')
a=re.findall(r'Total energy\s*:\s*([-+0-9.Ee]+)',out);assert len(a)==1,out
energy=float(a[0]);assert math.isfinite(energy),out
out=run('predict','predict.in','3.xsf');b=float(re.findall(r'Total energy\s*:\s*([-+0-9.Ee]+)',out)[0]);assert abs(energy-b)<1e-8
print('generate/train/predict passed; repeated energy:',energy)
