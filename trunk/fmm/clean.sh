find . \( -name '*.java' -o -name '*.class' -o -name '*.pem' -o -name '*_x10stub.c' -o -name 'javacore*' -o -name 'heapdump*' -o -name '*.stackdump' -o -name 'log*' -o -name 'TMP_*' \) -print0 | xargs -r0 rm -f

