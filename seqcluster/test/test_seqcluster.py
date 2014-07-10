from unittest import TestCase

import seqcluster


class TestJoke(TestCase):
    def test_prepare(self):
        s = seqcluster.prepare()
        self.assertTrue(isinstance(s, basestring))