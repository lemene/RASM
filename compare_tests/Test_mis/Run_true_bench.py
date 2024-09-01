
import os
import sys
work_dir = ""
# tool_ls = ["flye", "wtdbg2", "hifiasm", "smartdenovo"]
# data_type = "hifi"
tool_ls = ["flye", "wtdbg2", "smartdenovo"]
data_type = "ont"   # 
for tool in tool_ls:
    tool_dir = work_dir + "/" + tool
    if not os.path.isdir(tool_dir): os.makedirs(tool_dir)
    # 
    ref = "data" + "/" + tool + "_" + data_type + ".fa"
    
