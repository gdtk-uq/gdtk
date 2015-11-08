#!/usr/bin/env python
import os, sys
RUN_CMD = 'e4shared --job=connection-test --run'
of = open('test-results.txt', 'w')

faces = ["north", "south", "east", "west", "top", "bottom"]
nOrnt = 4
count = 0
for fA in faces:
    for fB in faces:
        for ornt in range(nOrnt):
            count += 1
            f = open('case.txt', 'w')
            f.write("%s\n" % fA)
            f.write("%s\n" % fB)
            f.write("%d\n" % ornt)
            f.close()
            # Prep stage
            val = os.system('e4shared --job=connection-test --prep')
            if val != 0:
                print "Problem preparing case: ", fA, fB, ornt
                sys.exit(1)
            # Run stage
            val = os.system(RUN_CMD)
            if val != 0:
                print "Problem running case: ", fA, fB, ornt
                sys.exit(1)
                # Post-process stage
            probe_str = ""
            for j in range(4):
                for k in range(4):
                    probe_str += "%f,%f,%f;" % (1.125, 0.125+0.25*j, 0.125+0.25*k)

            cmd = 'e4shared --job=connection-test --post --tindx-plot=last --probe="%s" --output-file=post.out' % probe_str[:-1]
            val = os.system(cmd)
            if val != 0:
                print "Problem post-processing case: ", fA, fB, ornt
                sys.exit(1)
            # Check if test passed
            cmd = 'diff ref.out post.out'
            val = os.system(cmd)
    
            of.write('CASE %03d: faceA=%s, faceB=%s, orientation=%d  \t' % 
                     (count, fA, fB, ornt))
            if val == 0:
                of.write('PASSED\n')
            else:
                of.write('FAILED\n')

of.close()



    
    


