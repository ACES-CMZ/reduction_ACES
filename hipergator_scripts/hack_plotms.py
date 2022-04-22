import sys
import os

apppath_string = '            app_path = __os.path.join( __os.path.abspath( __os.path.join(__os.path.dirname(__file__),"..") ), \'__bin__/casaplotms-x86_64.AppImage\')\n'
new_apppath_string = '            app_path = "/tmp/casaplotms/squashfs-root/AppRun"\n'

assert 'plotmstool.py' in sys.argv[1]
plotmstool = sys.argv[1]

with open(plotmstool, 'r') as fh:
    lines = fh.readlines()

newstr = """
            app_path = "/tmp/casaplotms/squashfs-root/AppRun"
            if not __os.path.exists(app_path):
                print(f"Did not find extracted path {app_path}")
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
else:
    print("Already modified.")
    ind = lines.index(new_apppath_string)
    print(f"Lines {ind}-{ind+4}: ")
    print(lines[ind:ind+4])

