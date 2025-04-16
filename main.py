"""
This main.py is new. Delete at a later time; I think.
2025-04-10.  Attempting to figure out and fix python importing.
    2025-04-16. Got import working, but the scheme is UGLY!
"""
import os

# from Python_Scripts.twobody3d import twobody3d
from Python_Test_Scripts.interplanetary_test import interplanetary_test

# the commented out code was used for testing imports
if __name__ == '__main__':
    current_directory = os.getcwd()
    print(f"Current working directory: {current_directory}")
    # twobody3d()
    interplanetary_test()
else:
    print("Other call method:")
