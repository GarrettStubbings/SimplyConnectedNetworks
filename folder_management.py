"""
Some stuff for dealing with folders
"""

import os
import time

def make_directory(folder):

    folders = [f for f in folder.split("/") if f != ""]

    if folder[0] == "/":
        path = "/"
    else:
        path = ""

    for f in folders:
        path += f + "/"
        if not os.path.isdir(path):
            print("Making:", path)
            os.mkdir(path)

def remove_directory(folder):
    folders = [f for f in folder.split("/") if f != ""]
    if folder[0] == "/":
        path = "/"
    else:
        path = ""

    for i in range(len(folders)):
        temp_path = path
        for f in folders[:len(folders)-i]:
            temp_path += f + "/"
        if temp_path[-8:] == "scratch/":
            print("Don't try to remove scratch.")
            break
        print("Removing:", temp_path)
        try:
            if os.path.isdir(temp_path):
                os.rmdir(temp_path)
            else:
                print(temp_path, "DNE?")
        except OSError:
            print("Non-Empty Directory:", temp_path)
            extra_files = os.listdir(temp_path)
            if len(extra_files) == 1 and extra_files[0] == ".DS_Store":
                print("Just the stupid DS_Store")
                os.remove(temp_path + ".DS_Store")
                os.rmdir(temp_path)
            else:
                print("Additional files, stopping.")
                break

