import unittest
import numpy as np
import scipy.sparse as sps


class TestNSETensor(unittest.TestCase):

    def setUp(self):
        self.NV = 100
        self.hmat = sps.rand(self.NV, self.NV*self.NV, format='csc')
        self.v = np.random.rand(self.NV, 1)

    def test_hmatvv(self):
        ''' evaluation of `h*vv` w/o forming `kron(v, v)` '''

        from conv_tensor_utils import eva_quadterm

        hvv = eva_quadterm(self.hmat, self.v)
        dhvv = self.hmat*np.kron(self.v, self.v)

        self.assertTrue(np.allclose(hvv, dhvv))

if __name__ == '__main__':
    unittest.main()
