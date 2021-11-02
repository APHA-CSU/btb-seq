import unittest
from subprocess import CompletedProcess
from unittest.mock import Mock, mock_open, patch

import validator

class ValidatorTests(unittest.TestCase):
    @patch("validator.open", new_callable=mock_open, read_data="mock-text")
    def test_validator(self, open_mock):
        """ A simple introductory unit test that checks we can run through the
            validator without errors
        """
        validator.subprocess.run = Mock(return_value=CompletedProcess([], 0))
        validator.subprocess.returncode = 0
        validator.open = mock_open()
        validator.os.makedirs = Mock()
        validator.os.path.isdir = Mock(side_effect=[False, True])
        validator.glob.glob = Mock(return_value=['./'])
        validator.analyse = Mock(return_value={})

        validator.performance_test('./', './', './reference.fasta')


if __name__ == '__main__':
    unittest.main()