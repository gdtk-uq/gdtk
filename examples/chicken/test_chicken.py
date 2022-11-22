"""
Integration tests for chicken, forward facing step edition

"""

import unittest
import chkn_post as chkn
import subprocess 
import os
from numpy import array, argmax

class TestFFS(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        os.chdir("forward-facing-step")
        cmd = subprocess.run('sed -i "s/config.max_time = 5.0e-3/config.max_time = 1.0e-3/" ffs.py', shell=True, stdout=subprocess.DEVNULL)

        result = subprocess.run("chkn-run", stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        result = result.stdout.decode('utf-8')
        version = result.splitlines()[1].split()[-1]

        result = subprocess.run("which chkn-run | xargs -I{} date -r {}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        date = result.stdout.decode('utf-8').strip()
        print("{} with chicken revision {}, installed {}".format(cls.__name__, version, date))

    def test_0_prep(self):
        cmd = subprocess.run('rm -rf ffs', shell=True, check=True, stdout=subprocess.DEVNULL)
        cmd = subprocess.run('chkn-prep -f ffs', shell=True, check=True, stdout=subprocess.DEVNULL)

    def test_1_run(self):
        cmd = subprocess.run('chkn-run --job=ffs', shell=True, check=True, stdout=subprocess.DEVNULL)

    def test_2_post(self):
        cmd = subprocess.run('chkn-post --job=ffs', shell=True, check=True, stdout=subprocess.DEVNULL)

    def test_3_shock_location(self):
        jobDir = "ffs"
        chkn.read_config(jobDir)
        flowblk000 = chkn.read_block_of_flow_data('ffs/flow/t0001/flow-0000-0000-0000.zip')
        nic = chkn.config["nics"][0]
        njc = chkn.config["njcs"][0]
        nkc = chkn.config["nkcs"][0]
        xblk = array(flowblk000["pos.x"]).reshape((nkc,njc,nic))
        pblk = array(flowblk000["p"]).reshape((nkc,njc,nic))
        p = pblk[0,0,:]
        x = xblk[0,0,:]

        # Search for the cell where the pressure increases above 1% of
        # the freestream. If this fails, shockidx will be 0.
        p0 = p[0]
        dp = 0.01*p0
        shockidx = argmax(p>p0+dp)
        self.assertNotEqual(shockidx,0)

        # We interpolate to find the position of the theoretical
        # normal shock pressure. This gives us an analogue shock 
        # position to test against.
        xp = x[shockidx-1]; pp = p[shockidx-1]
        xs = x[shockidx]; ps = p[shockidx]
        xss= x[shockidx+1]; pss= p[shockidx+1]

        pn = 1047025.00
        xn = (xs-xp)/(ps-pp)*(pn-pp) + xp

        # one tenth of the cell size seems like a pretty tight tolerance
        tol = (xs-xp)*0.1 

        # Arbitrary, from commit 31cf2ca84f1e0452be73c015010539c930a8fa2e
        shock_position_target = 0.41369441539674984
        self.assertEqual(True, abs(xn-shock_position_target)<tol)
        
        # Only clean up with the tests have worked correctly
        cmd = subprocess.run('rm -rf ffs', shell=True, check=True, stdout=subprocess.DEVNULL)

    @classmethod
    def tearDownClass(cls):
        cmd = subprocess.run('sed -i "s/config.max_time = 1.0e-3/config.max_time = 5.0e-3/" ffs.py', shell=True, stdout=subprocess.DEVNULL)
        os.chdir("../")

if __name__=='__main__':
    unittest.main()
