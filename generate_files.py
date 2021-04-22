import sys
from shutil import copyfile



def add_info(i, file):
    new_content = "// "+"p"+str(i+1)+".cpp"+"\n"
    f = open(file, 'r+')
    lines = f.readlines()  # read old content
    f.seek(0)  # go back to the beginning of the file
    f.write(new_content)  # write new content at the beginning
    for line in lines:  # write old content after new
        f.write(line)
    f.close()

n = int(sys.argv[1])
for i in range(n):
    file_name = "/Users/ashishsingh/Documents/cpp/codeforces/"+"p"+str(i+1)+".cpp"
    src_file = "/Users/ashishsingh/Documents/cpp/codeforces/template.cpp"
    copyfile(src_file, file_name)
    add_info(i, file_name)