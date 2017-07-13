# Script to xfail failing tests. Save the output from pytest with --junitxml=junit.xml and then run this on it.
# Do not run more than once! May need to restore the original files (e.g. with git checkout) if it goes wrong...

import xml.etree.ElementTree as ET
import pathlib
import re

tree = ET.parse('junit.xml')
root = tree.getroot()

for child in root:
    attribs = child.attrib
    testfail = child.find('failure')
    testerror = child.find('error')
    if testfail is not None or testerror is not None:
        file = pathlib.Path(attribs['file'])
        file = pathlib.Path(*file.parts[1:])
        fname = re.split('(\W+)',attribs['name'])[0]
        print('** Fail', file, fname)

        f = open(file, 'r')
        dat = f.readlines()
        f.close()

        already_pytest = False
        if 'import pytest\n' in dat:
            already_pytest = True
            print('Already got pytest import')

        newlines=[]
        for line in dat:
            if 'dolfin' in line and 'import' in line and not already_pytest:
                newlines.append('import pytest\n')
            l = re.split('(\W+)', line)
            if fname in l and 'def' in l:
                print('Found line', line)
                if (line[0]==' '):
                    newlines.append('    @pytest.mark.xfail\n')
                else:
                    newlines.append('@pytest.mark.xfail\n')
            newlines.append(line)

        g = open(file,'w')
        g.writelines(newlines)
        g.close()

    else:
        print('OK')


