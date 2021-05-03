import sys
from shutil import copyfile



def add_info(name, file):
    new_content = "// " + name + ".cpp" + "\n"
    f = open(file, 'r+')
    lines = f.readlines()  # read old content
    f.seek(0)  # go back to the beginning of the file
    f.write(new_content)  # write new content at the beginning
    for line in lines:  # write old content after new
        f.write(line)
    f.close()

n = int(sys.argv[1])
for i in range(n):
    name = "p" + str(i+1)
    file_name = "../CurrentContest/" + name + ".cpp"
    src_file = "../src/template.cpp"
    copyfile(src_file, file_name)
    add_info(name, file_name)