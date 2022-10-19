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

n = str(sys.argv[1])
subfolder = "PS_Codeforces"
if(len(sys.argv) > 2):
    subfolder = str(sys.argv[2])
template = "template_basic"
# print(len(sys.argv))
if(len(sys.argv) > 3):
    template = str(sys.argv[3])
file_name = "../ProblemSet/" + subfolder + "/" + n + ".cpp"
src_file = "../src/" + template + ".cpp"
copyfile(src_file, file_name)
add_info(n, file_name)