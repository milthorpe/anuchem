#!/bin/bash
main() {
  TESTLIB=~/x10-17/x10.tests/examples/x10lib:../src
  for file in `find . -name *.x10`
  do
    x10c -v -sourcepath $TESTLIB $file
    local CLASSPATH="`sed -ne 's|^\s*//\s*CLASSPATH*\:\s*\(.*\)|\1|p' "$file"`"
    local test="${file%.x10}"
    removeDotSlash "$test" 
    changeSlashToDot "$RESULT"
    local CLASSNAME="$RESULT"
    x10 -t -v -mx 128M -classpath \"$testDir\" $CLASSPATH $CLASSNAME
  done
}

removeDotSlash() {
# convert ./a/b/./c to a/b/c
  local str=$1
  #str=${str//\.\//}
  str=`echo "$str" | sed -e 's/\.\///g'`

  RESULT=$str
}

changeSlashToDot() {
# convert a/b/c to a.b.c, but ./foo becomes foo
  local str=$1
  str=${str#\.\/}
  #str=${str//\//\.}
  str=`echo "$str" | sed -e 's/\//\./g'`
  RESULT="$str"
}

CURRENT_DIR=`pwd`
cd ~/testx10/fmm/test
main
cd $CURRENT_DIR
exit 0

