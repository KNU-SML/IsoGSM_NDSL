#!/bin/csh -f
ps aux | egrep -i "kana|kei|jritchie|ldehaan|elena" | egrep 'atlas-|compas-' | sed '/grep/d' | awk '{print $1,$12}' | uniq -c

