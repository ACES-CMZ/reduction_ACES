import sys
import os
import tempfile
tmpdir = tempfile.gettempdir()

print(f"Using casaplotms from directory {tmpdir}")
plotmspath = tmpdir+'/casaplotms/squashfs-root/AppRun'
print(f"Does the path {plotmspath} exist?: {os.path.exists(plotmspath)}")

apppath_string = '            app_path = __os.path.join( __os.path.abspath( __os.path.join(__os.path.dirname(__file__),"..") ), \'__bin__/casaplotms-x86_64.AppImage\')\n'
new_apppath_string = '            app_path = "'+tmpdir+'/casaplotms/squashfs-root/AppRun"\n'

assert 'plotmstool.py' in sys.argv[1]
plotmstool = sys.argv[1]

with open(plotmstool, 'r') as fh:
    lines = fh.readlines()

if len(lines) < 900:
    print("plotms is corrupted!  re-load it!")
    sys.exit(99)

newstr = f"""
            {new_apppath_string}
            if not __os.path.exists(app_path):
                print(f"Did not find extracted path {{app_path}}")
                app_path = __os.path.join( __os.path.abspath( __os.path.join(__os.path.dirname(__file__),"..") ), '__bin__/casaplotms-x86_64.AppImage')
"""

if apppath_string in lines:
    linenum = lines.index(apppath_string)

    with open(plotmstool, 'w') as fh:
        for line in lines[:linenum]:
            fh.write(line)
        fh.write(newstr)
        for line in lines[linenum+1:]:
            fh.write(line)
    print(f"Hacked line {linenum} of {plotmstool}")
elif new_apppath_string not in lines:

    # find the index of the outdated app_path
    for ind, line in enumerate(lines):
        if 'app_path = "' in line:
            print("app_path was previously set to:")
            print(line)
            break

    linenum = ind

    with open(plotmstool, 'w') as fh:
        for line in lines[:linenum]:
            fh.write(line)
        fh.write(new_apppath_string)
        for line in lines[linenum+1:]:
            fh.write(line)
    print(f"Hacked line {linenum} of {plotmstool}")

else:
    print("Already modified.")
    linenum = lines.index(new_apppath_string)

print(f"Lines {linenum}-{linenum+4}: ")
print("".join(lines[linenum:linenum+4]))

