import os

files = os.listdir()
for file in files:
    try:
        with open(file, "r") as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        continue
    except IsADirectoryError:
        continue
    else:
        if lines == [] or lines == ["\n"]:
            print(f"Removing empty file: {file}")
            os.remove(file)
