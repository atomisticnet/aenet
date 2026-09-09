# Experimental AENET packaging; MPL-2.0, see src/license-header.txt.
import pathlib, platform, re, shutil, subprocess, sys
root=pathlib.Path(sys.argv[1]).resolve();mac=platform.system()=='Darwin'
def out(*args):return subprocess.check_output(args,text=True)
queue=[p for d in ['bin','tools','lib'] for p in (root/d).iterdir() if p.is_file() and not p.name.endswith('.a')]
seen=set()
while queue:
 p=queue.pop()
 if p in seen:continue
 seen.add(p)
 deps=([l.strip().split(' (')[0] for l in out('otool','-L',str(p)).splitlines()[1:]] if mac else re.findall(r'=>\s+(/\S+)',out('ldd',str(p))))
 for dep in deps:
  if mac:
   if not dep.startswith('/') or dep.startswith(('/usr/lib/','/System/')):continue
  elif pathlib.Path(dep).name.startswith(('libc.so','libm.so','libpthread.so','libdl.so','librt.so')):continue
  dest=root/'lib'/pathlib.Path(dep).name
  if not dest.exists():shutil.copy2(dep,dest);queue.append(dest)
  if mac:
   replacement=('@loader_path/' if p.parent.name=='lib' else '@loader_path/../lib/')+dest.name
   subprocess.run(['install_name_tool','-change',dep,replacement,str(p)],check=True)
 if mac:
  if p.suffix=='.dylib':subprocess.run(['install_name_tool','-id','@rpath/'+p.name,str(p)],check=True)
 else:subprocess.run(['patchelf','--set-rpath','$ORIGIN' if p.parent.name=='lib' else '$ORIGIN/../lib',str(p)],check=True)
for p in sorted(seen):
 if mac:subprocess.run(['codesign','--force','--sign','-',str(p)],check=True)
 print(p.relative_to(root));print(out('otool','-L',str(p)) if mac else out('ldd',str(p)))
