import os
import shutil


for i in xrange(101, 170):
    src = "sim003" + str(i) + ".out"
    dst = "sim003" + str(i)
    shutil.move(src, dst)



for i in xrange(101, 170):
    main_dir = "sim003" + str(i)
    for j in xrange(0, 7):
       src = main_dir + "_0" + str(j)
       dst = main_dir
       shutil.move(src, dst)
