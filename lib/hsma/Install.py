#!/usr/bin/env python
#coding:utf-8

"""
Install.py tool to download, unpack, build, and link to the HSMA library
used to automate the steps described in the README file in this dir
"""

from __future__ import print_function
import sys, os, subprocess, shutil, zipfile
from argparse import ArgumentParser

sys.path.append('..')
from install_helpers import fullpath, geturl, get_cpus, checkmd5sum

parser = ArgumentParser(prog='Install.py', 
                        description="LAMMPS library build wrapper script")

# settings

version = "1.0.0" 
url = "https://github.com/LiangJiuyang/FMM3D/blob/master/LibHSMA/HSMA%s.zip" % (version)

# extra help message

HELP = """
Syntax from src dir: make lib-hsma args="-b"
                 or: make lib-hsma args="-p /usr/local/hsma"
Syntax from lib dir: python Install.py -b
                 or: python Install.py -p /usr/local/hsma

Example:

make lib-hsma args="-b"   # download/build in lib/hsma/hsma
make lib-hsma args="-p $HOME/hsma" # use existing HSMA installation in $HOME
"""

# parse and process arguments

pgroup = parser.add_mutually_exclusive_group()
pgroup.add_argument("-b", "--build", action="store_true",
                    help="download and build the HSMA library") 
pgroup.add_argument("-p", "--path",
                    help="specify folder of existing HSMA installation")
parser.add_argument("-v", "--version", default=version,
                    help="set version of HSMA to download and build (default: %s)" % version)

args = parser.parse_args()

# print help message and exit, if neither build nor path options are given
if not args.build and not args.path:
  parser.print_help()
  sys.exit(HELP)

buildflag = args.build
pathflag = args.path is not None
version = args.version

homepath = fullpath(".")
hsmapath = os.path.join(homepath, "HSMA%s" % version)

#print(scafacospath)

#print(scafacospath)

if pathflag:
  hsmapath = args.path 
  if not os.path.isdir(os.path.join(hsmapath, "include")):
    sys.exit("HSMA include path for %s does not exist" % hsmapath) 
  if (not os.path.isdir(os.path.join(hsma, "lib64"))) \
     and (not os.path.isdir(os.path.join(hsmapath, "lib"))):
    sys.exit("HSMA lib path for %s does not exist" % hsmapath) 
  hsmapath = fullpath(hsmapath)

#print(fullpath(scafacospath))

# download and unpack tarball 

if buildflag: 
  print("Downloading HSMA ...") 
  
  print(hsmapath)
  print(url)
  print(homepath)

  geturl(url, "%s/HSMA%s.zip" % (homepath, version))

  print("Unpacking HSMA tarball ...")

  if os.path.exists(hsmapath): 
    shutil.rmtree(hsmapath) 

  zipname = os.path.join(homepath, "%s.zip" % hsmapath) 
  print(zipname)

  if zipfile.is_zipfile(zipname): 
    tgz = zipfile.open(zipname) 
    tgz.extractall(path=homepath) 
    os.remove(zipname)
  else:
    cmd = 'unzip HSMA%s.zip' % (version)
    txt = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True) 
    #sys.exit("File %s is not a supported archive" % zipname) 


print("Creating links to HSMA include and lib files")

if os.path.isfile("includelink") or os.path.islink("includelink"):
  os.remove("includelink")

if os.path.isfile("liblink") or os.path.islink("liblink"):
  os.remove("liblink")

if buildflag:
  os.symlink(os.path.join(homepath, 'build', 'include'), 'includelink')
  os.symlink(os.path.join(homepath, 'build', 'lib'), 'liblink')
else:
  os.symlink(os.path.join(hsmapath, 'include'), 'includelink')
  if os.path.isdir(os.path.join(hsmapath, "lib64")):
    os.symlink(os.path.join(hsmapath, 'lib64'), 'liblink')
  else:
    os.symlink(os.path.join(hsmapath, 'lib'), 'liblink')