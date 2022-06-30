import unittest
import subprocess

class GetSubmissionNoTests(unittest.TestCase):
    
    def test_get_submission_no_test(self):
        test_input = ["AFxx-12-34567-89",
                      "ATxx-12-34567-89",
                      "AFx-12-34567-89",
                      "Ax-12-34567-89",
                      "AF-12-34567-89",
                      "AFx12-34567-89",
                      "HI-12-34567-89",
                      "12-34567-89-1L",
                      "12-34567-89-L1",
                      "A-12-34567-89",
                      "12-34567-89-1",
                      "12-34567-89-L",
                      "12-34567-89",
                      "AFxx-12-3456-89",
                      "ATxx-12-3456-89",
                      "AFx-12-3456-89",
                      "Ax-12-3456-89",
                      "AF-12-3456-89",
                      "AFx12-3456-89",
                      "HI-12-3456-89",
                      "12-3456-89-1L",
                      "12-3456-89-L1",
                      "A-12-3456-89",
                      "12-3456-89-1",
                      "12-3456-89-L",
                      "12-3456-89",
                      "12345678",
                      "ABCDEFGH",
                      ""]
        test_output = ["12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-34567-89",
                       "12-3456-89",
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12-3456-89", 
                       "12345678", 
                       "ABCDEFGH", 
                       ""]
        fail = False 
        i = 0
        for n, (input, output) in enumerate(zip(test_input, test_output)):
            try:
                ps = subprocess.run(["bash", "./bin/get_submission_no.bash", input], capture_output=True)
                self.assertEqual(ps.stdout.decode().strip('\n'), output)
            except AssertionError as e:
                i += 1
                fail = True
                print(f"Test failure #{n}: ", e)
        if fail: 
            print(f"{i} test failures")
            raise AssertionError

if __name__ == "__main__":
    unittest.main()
